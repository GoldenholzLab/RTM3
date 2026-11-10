#!/usr/bin/env python3
"""Export the executed appendix to self-contained HTML and native-math DOCX.

No simulation, office renderer, journal upload, or repository write is invoked.
Requires pandoc plus the Python reporting/document dependencies.
"""
import argparse
import base64
import copy
import hashlib
import json
import platform
import re
import subprocess
from pathlib import Path
from zipfile import ZipFile

import nbformat
from bs4 import BeautifulSoup
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK, WD_TAB_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor
from lxml import etree

ROOT = Path(__file__).resolve().parents[1]
BASE_COMMIT = "59b47a5564342d4b2695307a341fbd0ff1d90091"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
      "m": "http://schemas.openxmlformats.org/officeDocument/2006/math"}
FIGURES = [
    ("appendix_fixed_vs_reselected.png", "Fixed versus reselected cohorts",
     "The lower-left panel uses 0 to 40 percent with ticks every 10 percentage points. Selection uses 100,000 potential participants; observed RTM uses all selected participants; MPC uses selected participants with positive measured baseline. The fixed cohort is a hypothetical diagnostic, as qualified above."),
    ("appendix_original_25_day.png", "Sensitivity analysis of the 25-day eligibility rule",
     "Only the interval criterion differs between the two series. Monthly counts and the reference retain the primary definitions. Selection and outcome denominators are those in the complete tables above."),
    ("appendix_correction_alternatives.png", "False-alarm correction and interval alternatives",
     "Sensitivity is 100%. Strategies share detector realizations but differ in interval definition or count correction and enrollment. Raw outcomes use the observation-scale reference including expected alarms. Oracle removal requires known simulated alarm identities."),
]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def table_rows(html):
    table = BeautifulSoup(html, "html.parser").find("table")
    head_rows = table.select("thead tr")
    body = [[c.get_text().strip() for c in tr.find_all(["th", "td"], recursive=False)]
            for tr in table.select("tbody tr")]
    head = [c.get_text().strip() for c in head_rows[0].find_all(["th", "td"], recursive=False)]
    # Pandas emits the named index in a second header row for the threshold table.
    if len(head_rows) == 2 and not head[0]:
        head[0] = head_rows[1].find(["th", "td"]).get_text().strip()
    return head, body


def markdown_table(head, body):
    def esc(x):
        return x.replace("|", "\\|").replace("\n", " ")
    return "\n".join(["| " + " | ".join(map(esc, head)) + " |",
                      "| " + " | ".join(["---"] * len(head)) + " |"] +
                     ["| " + " | ".join(map(esc, row)) + " |" for row in body])


def assembled_source(nb, out):
    pieces, code_blocks, inventory = [], [], []
    section, group, table_number, figure_number = "Design", "", 0, 0
    for ci, cell in enumerate(nb.cells):
        if cell.cell_type == "markdown":
            pieces.append(cell.source)
            for line in cell.source.splitlines():
                if line.startswith("## "):
                    section, group = line[3:], ""
            inventory.append({"cell": ci, "kind": "narrative", "sha256": sha(cell.source.encode())})
            continue
        code_blocks.append((cell.execution_count, cell.source))
        for oi, output in enumerate(cell.outputs):
            data = output.get("data", {})
            item = {"cell": ci, "output": oi}
            if "image/png" in data:
                name, title, note = FIGURES[figure_number]
                image_bytes = base64.b64decode(data["image/png"])
                (ROOT / "results" / name).write_bytes(image_bytes)
                image_path = out / name
                image_path.write_bytes(image_bytes)
                figure_number += 1
                pieces.append(f"![Figure {figure_number}. {title}.]({name}){{width=9.4in}}\n\n{note}")
                item.update(kind="figure", number=figure_number, filename=name, sha256=sha(image_bytes))
            elif "text/html" in data:
                head, rows = table_rows(data["text/html"])
                table_number += 1
                title = section + (" — " + group if group else "")
                if table_number == 1:
                    title = "Recorded simulation design"
                elif table_number == 2:
                    title = "Recorded simulation software"
                elif "Retained n" in head:
                    title = "Fixed membership and current detector eligibility"
                elif "each_baseline_month_threshold_raw" in head:
                    title = "Raw additive thresholds and eligibility mask agreement"
                elif "Maximum absolute difference, pp" in head:
                    title = "Largest absolute differences after retaining the original 25-day rule"
                elif "Selection difference, pp" in head:
                    title = "Original 25-day rule minus primary analysis"
                pieces.append(f"Table {table_number}. {title}\n\n" + markdown_table(head, rows))
                item.update(kind="table", number=table_number, title=title, headers=head, rows=rows,
                            body_row_count=len(rows))
            elif "text/markdown" in data:
                value = data["text/markdown"]
                pieces.append(value)
                if value.startswith("### "):
                    group = value[4:].strip()
                item.update(kind="markdown_output", sha256=sha(value.encode()))
            elif "text/plain" in data:
                value = data["text/plain"]
                pieces.append("```text\n" + value + "\n```")
                item.update(kind="plain_output", text=value)
            elif output.output_type == "stream":
                pieces.append("```text\n" + output.text + "\n```")
                item.update(kind="stream", text=output.text)
            else:
                raise ValueError(f"Unaccounted output {ci}:{oi}")
            inventory.append(item)
    pieces.append("## Reporting code and execution\n\n"
                  f"Source repository: <https://github.com/GoldenholzLab/RTM3/tree/{BASE_COMMIT}>. "
                  "The following are all code cells from this executed appendix, in execution order. "
                  "The notebook delegates reporting to rtm3_reporting.py. The remote reporting run reads the released CSV and JSON files "
                  "and the saved results of the new paired analysis, checks source and artifact hashes, and regenerates figures. "
                  "Notebook execution does not rerun the simulation. The separate, authorized participant-level replay is described "
                  "in the paired-analysis methods and provenance sections. "
                  "The full run metadata above retains the original simulation provenance. "
                  "The export audit separately records reporting software and the native-equation conversion chain.\n\n"
                  "This remote host disallows TCP and IPC socket listeners. All code cells were therefore executed top to bottom "
                  "in a fresh in-process IPython shell with rich display capture and the inline plotting backend. "
                  "The saved notebook retains execution counts and all outputs.\n\n"
                  "To rebuild the notebook and exports from the repository root:\n\n"
                  "```bash\npython scripts/build_appendix_notebook.py --execute --kernel-transport inprocess\n"
                  "python scripts/export_appendix.py\n```\n\n"
                  "Standard Jupyter execution remains available with the builder's default TCP kernel transport on hosts that allow it.")
    for number, source in code_blocks:
        pieces.append(f"### Code cell {number}\n\n```python\n{source}\n```")
    return "\n\n".join(pieces) + "\n", inventory


def field(paragraph, instruction, cached):
    for kind, value in [("begin", None), ("instr", instruction), ("separate", None), ("text", str(cached)), ("end", None)]:
        run = paragraph.add_run()
        if kind == "instr":
            el = OxmlElement("w:instrText"); el.set(qn("xml:space"), "preserve"); el.text = value
        elif kind == "text":
            run.text = value
            continue
        else:
            el = OxmlElement("w:fldChar"); el.set(qn("w:fldCharType"), kind)
        run._r.append(el)


def style_document(path):
    doc = Document(path)
    for section in doc.sections:
        section.orientation = WD_ORIENT.LANDSCAPE
        section.page_width, section.page_height = Inches(11), Inches(8.5)
        section.top_margin = section.bottom_margin = Inches(.65)
        section.left_margin = section.right_margin = Inches(.7)
        section.header_distance = section.footer_distance = Inches(.28)
    styles = {s.name: s for s in doc.styles}
    for style in doc.styles:
        if style.type == 1 or style.type == 2:
            style.font.name = "Times New Roman"
            if style.type == 1: style.font.size = Pt(11.5)
            style.font.color.rgb = RGBColor(0, 0, 0)
    normal = styles["Normal"]
    normal.font.size = Pt(11.5)
    normal.paragraph_format.line_spacing = 1.12
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.widow_control = True
    for name, size in [("Title", 18), ("Heading 1", 15), ("Heading 2", 13), ("Heading 3", 11.5)]:
        style = styles[name]; style.font.size = Pt(size); style.font.bold = True
        style.paragraph_format.space_before = Pt(12); style.paragraph_format.space_after = Pt(6)
        style.paragraph_format.keep_with_next = True
    for name in ["Source Code", "Verbatim Char"]:
        if name in styles:
            styles[name].font.name = "Consolas"; styles[name].font.size = Pt(9.5)
    if "Image Caption" in styles:
        styles["Image Caption"].font.size = Pt(11)
    equation_number = 0
    for p in doc.paragraphs:
        if p.style.name == "Heading 1" and not any(x.style.name == "Title" for x in doc.paragraphs[:1]):
            p.style = styles["Title"]
        if p.text.startswith("Table "):
            p.paragraph_format.keep_with_next = True
            p.paragraph_format.space_before = Pt(10)
            p.paragraph_format.space_after = Pt(6)
            for r in p.runs: r.bold = True; r.font.size = Pt(11)
        if p._p.xpath(".//w:drawing"):
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.keep_with_next = True
            p.paragraph_format.space_before = Pt(10)
        # Pandoc uses m:oMathPara for displays. Keep native OMML and add live SEQ fields.
        displays = p._p.findall(qn("m:oMathPara"))
        if displays:
            equation_number += 1
            maths = [copy.deepcopy(m) for d in displays for m in d.findall(qn("m:oMath"))]
            for d in displays: p._p.remove(d)
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            p.paragraph_format.tab_stops.add_tab_stop(Inches(4.65), WD_TAB_ALIGNMENT.CENTER)
            p.paragraph_format.tab_stops.add_tab_stop(Inches(9.55), WD_TAB_ALIGNMENT.RIGHT)
            p.paragraph_format.keep_together = True
            p.paragraph_format.space_before = Pt(8); p.paragraph_format.space_after = Pt(8)
            bookmark = OxmlElement("w:bookmarkStart"); bookmark.set(qn("w:id"), str(5000+equation_number)); bookmark.set(qn("w:name"), f"Equation{equation_number}")
            p._p.append(bookmark)
            p.add_run("\t")
            for math in maths: p._p.append(math)
            p.add_run("\t("); field(p, " SEQ Equation \\* ARABIC ", equation_number); p.add_run(")")
            end = OxmlElement("w:bookmarkEnd"); end.set(qn("w:id"), str(5000+equation_number)); p._p.append(end)
    # Pandoc can size images for its initial portrait page before this exporter
    # applies landscape sections. Set the final width without changing PNG bytes.
    for shape in doc.inline_shapes:
        original_width, original_height = shape.width, shape.height
        final_width = Inches(9.4)
        shape.width = final_width
        shape.height = round(original_height * final_width / original_width)
    for table in doc.tables:
        table.autofit = False
        n = len(table.columns)
        if n == 8: widths = [1.42, .90, .85, 1.18, .83, 1.60, 1.10, .72]
        elif n == 9: widths = [1.25, .70, .65, 1.05, 1.25, .80, .80, 1.55, 1.55]
        elif n == 5: widths = [1.35, 1.00, 2.25, 2.25, 1.75]
        elif n == 2: widths = [4.5, 5.1]
        else: widths = [9.6/n]*n
        scale = 9.6/sum(widths); widths = [w*scale for w in widths]
        for col, width in zip(table.columns, widths): col.width = Inches(width)
        pr = table._tbl.tblPr
        borders = OxmlElement("w:tblBorders")
        for side in ["top", "left", "bottom", "right", "insideH", "insideV"]:
            el = OxmlElement("w:"+side); el.set(qn("w:val"), "single"); el.set(qn("w:sz"), "4"); el.set(qn("w:color"), "D9D9D9"); borders.append(el)
        for old in pr.findall(qn("w:tblBorders")): pr.remove(old)
        pr.append(borders)
        for ri, row in enumerate(table.rows):
            if ri == 0:
                props = row._tr.get_or_add_trPr()
                if props.find(qn("w:tblHeader")) is None:
                    repeat = OxmlElement("w:tblHeader"); props.append(repeat)
            props = row._tr.get_or_add_trPr()
            if props.find(qn("w:cantSplit")) is None:
                cant = OxmlElement("w:cantSplit"); props.append(cant)
            for ci, cell in enumerate(row.cells):
                cell.width = Inches(widths[ci]); cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
                cp = cell._tc.get_or_add_tcPr()
                margins = OxmlElement("w:tcMar")
                for side, val in [("top",85),("bottom",85),("left",90),("right",90)]:
                    el=OxmlElement("w:"+side);el.set(qn("w:w"),str(val));el.set(qn("w:type"),"dxa");margins.append(el)
                for old in cp.findall(qn("w:tcMar")): cp.remove(old)
                cp.append(margins)
                if ri == 0:
                    shade=OxmlElement("w:shd");shade.set(qn("w:fill"),"E7E7E7");cp.append(shade)
                for p in cell.paragraphs:
                    p.paragraph_format.space_before = Pt(0); p.paragraph_format.space_after = Pt(0)
                    p.paragraph_format.line_spacing = 1.05; p.paragraph_format.keep_with_next = ri == 0
                    p.alignment = WD_ALIGN_PARAGRAPH.LEFT if ci == 0 or n == 2 else WD_ALIGN_PARAGRAPH.CENTER
                    for run in p.runs:
                        run.font.size = Pt(10.5); run.bold = ri == 0; run.font.color.rgb = RGBColor(0,0,0)
                        # Permit Word to wrap long metadata/header identifiers at underscores.
                        if "_" in run.text: run.text = run.text.replace("_", "_\u200b")
    section=doc.sections[0]
    header=section.header.paragraphs[0]; header.text="RTM3 scientific appendix"
    header.runs[0].font.size=Pt(9)
    footer=section.footer.paragraphs[0];footer.alignment=WD_ALIGN_PARAGRAPH.RIGHT
    footer.add_run("Pending authoritative Microsoft Word verification  |  ").font.size=Pt(9)
    field(footer," PAGE ",1)
    settings=doc.settings.element
    for old in settings.findall(qn("w:updateFields")): settings.remove(old)
    update=OxmlElement("w:updateFields");update.set(qn("w:val"),"true");settings.append(update)
    doc.core_properties.title="RTM3 scientific appendix"
    doc.core_properties.subject="Detector error eligibility and observed regression to the mean"
    doc.core_properties.author=""; doc.core_properties.last_modified_by=""
    doc.core_properties.comments="Generated remotely from released aggregates; pending authoritative Microsoft Word verification."
    doc.save(path)
    return equation_number


def math_audit(ast, docx, out):
    equations=[]
    def walk(node):
        if isinstance(node,dict):
            if node.get("t")=="Math":
                equations.append({"kind":node["c"][0]["t"],"latex":node["c"][1]})
            for val in node.values(): walk(val)
        elif isinstance(node,list):
            for val in node:walk(val)
    walk(ast)
    with ZipFile(docx) as z: root=etree.fromstring(z.read("word/document.xml"))
    maths=root.xpath("//m:oMath",namespaces=NS)
    assert len(maths)==len(equations),(len(maths),len(equations))
    display_no=0
    for i,(source,omml) in enumerate(zip(equations,maths),1):
        source["math_index"]=i
        source["omml_file"]=f"math_{i:03}.xml"
        if source["kind"]=="DisplayMath": display_no+=1;source["live_equation_number"]=display_no
        (out/source["omml_file"]).write_bytes(etree.tostring(omml,encoding="utf-8",xml_declaration=True,pretty_print=True))
        source["omml_text"]="".join(omml.xpath(".//m:t/text()",namespaces=NS))
    fields=root.xpath("//w:instrText/text()",namespaces=NS)
    assert sum("SEQ Equation" in f for f in fields)==display_no
    (out/"equations.json").write_text(json.dumps(equations,indent=2,ensure_ascii=False))
    (out/"equations.tex").write_text("% Exact LaTeX from the assembled appendix source\n"+"\n\n".join("% Equation "+str(e["live_equation_number"])+"\n\\[\n"+e["latex"]+"\n\\]" for e in equations if e["kind"]=="DisplayMath"))
    return {"inline_and_display_math_count":len(equations),"display_equations":display_no,
            "live_SEQ_equation_fields":display_no,"converter":"pandoc LaTeX reader / texmath to OMML writer",
            "chain":"RTM3_appendix_source.md -> appendix_ast.json -> DOCX word/document.xml; equations.json maps every LaTeX expression to math_NNN.xml",
            "authoritative_Word_equation_visual_verification":"PENDING"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,default=ROOT/"results"/"appendix_export")
    parser.add_argument("--html-only",action="store_true",
                        help="Export notebook/HTML and source audits only; leave every existing DOCX untouched.")
    args=parser.parse_args();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    audit_dir=out/"native_equations";audit_dir.mkdir(exist_ok=True)
    nb=nbformat.read(ROOT/"for_appendix_RTM3.ipynb",as_version=4);nbformat.validate(nb)
    counts=[c.execution_count for c in nb.cells if c.cell_type=="code"]
    assert counts==list(range(1,len(counts)+1)),counts
    assert all(o.output_type!="error" for c in nb.cells for o in c.get("outputs",[]))
    source,inventory=assembled_source(nb,out)
    source_file=out/"RTM3_appendix_source.md";source_file.write_text(source)
    # Keep the LaTeX source and Pandoc AST, not an uninspectable equation image.
    ast_result=subprocess.run(["pandoc",str(source_file),"-f","markdown","-t","json"],capture_output=True,text=True,check=True)
    ast=json.loads(ast_result.stdout);(audit_dir/"appendix_ast.json").write_text(ast_result.stdout)
    docx=out/"RTM3_appendix.docx"
    if not args.html_only:
        result=subprocess.run(["pandoc",str(source_file),"--resource-path",str(out),"-f","markdown","-t","docx","-o",str(docx)],capture_output=True,text=True,check=True)
        (audit_dir/"conversion.log").write_text(result.stderr)
        assert not result.stderr.strip(),result.stderr
        displays=style_document(docx)
    css=out/"appendix.css"
    css.write_text("body{font-family:Georgia,serif;max-width:1280px;margin:2em auto;padding:0 2em;color:#111;font-size:17px;line-height:1.5}h1,h2,h3{color:#000;line-height:1.25}h2{margin-top:2em}table{border-collapse:collapse;width:100%;font-size:14px;margin:1em 0 2em}th,td{padding:9px 8px;border:1px solid #d9d9d9;text-align:center;overflow-wrap:anywhere}th{background:#e7e7e7}td:first-child,th:first-child{text-align:left}thead{display:table-header-group}img{max-width:100%;height:auto}figure{margin:2em 0}figcaption{font-size:16px;font-weight:bold}pre{background:#f7f7f7;padding:1em;white-space:pre-wrap;overflow-wrap:anywhere;font-size:13px}math[display=block]{margin:1em auto;overflow-x:auto}a{color:#123d67}@media print{body{max-width:none;margin:0}table{font-size:10pt}tr{break-inside:avoid}h2,h3,figcaption{break-after:avoid}}")
    html=ROOT/"results"/"notebooks"/"for_appendix_RTM3.html"
    subprocess.run(["pandoc",str(source_file),"--resource-path",str(out),"-f","markdown","-t","html5","--standalone","--embed-resources","--mathml","--css",str(css),"--metadata","title=RTM3 scientific appendix","--template",str(write_html_template(out)),"-o",str(html)],check=True)
    (out/html.name).write_bytes(html.read_bytes())
    (out/"for_appendix_RTM3.ipynb").write_bytes((ROOT/"for_appendix_RTM3.ipynb").read_bytes())
    if not args.html_only:
        math_result=math_audit(ast,docx,audit_dir)
        assert displays==math_result["display_equations"]
    else:
        math_result={"status":"Existing Word equations untouched; HTML-only export. Await current user v1 DOCX before direct OOXML insertion.",
                     "source_chain":"RTM3_appendix_source.md -> native_equations/appendix_ast.json -> HTML MathML"}
    import importlib.metadata as md
    packages=["numpy","pandas","scipy","matplotlib","nbformat","nbclient","ipython","python-docx","beautifulsoup4","lxml","jinja2","matplotlib-inline","ipykernel","jupyter-client"]
    audit={"base_repository_commit":BASE_COMMIT,"reporting_platform":platform.platform(),"reporting_python":platform.python_version(),"reporting_software":{p:md.version(p) for p in packages},
           "pandoc":subprocess.check_output(["pandoc","--version"],text=True).splitlines()[0],
           "executed_code_cells":counts,"simulation_rerun_in_notebook_or_export":False,
           "separate_paired_replay":"Recorded in results/paired_type3/paired_metadata.json",
           "recorded_simulation_metadata_modified":False,
           "tables":len([i for i in inventory if i.get("kind")=="table"]),"figures":3,
           "equations":math_result,"word_status":"PENDING authoritative Microsoft Word verification; no Office renderer invoked",
           "html_math":"Native MathML; all figure images and CSS embedded; no external rendering scripts"}
    (out/"reporting_audit.json").write_text(json.dumps(audit,indent=2))
    (out/"content_inventory.json").write_text(json.dumps(inventory,indent=2,ensure_ascii=False))
    print(json.dumps(audit,indent=2))


def write_html_template(out):
    template=out/"appendix_template.html"
    template.write_text('<!DOCTYPE html>\n<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>$pagetitle$</title>$for(css)$<link rel="stylesheet" href="$css$">$endfor$</head><body><main>$body$</main></body></html>')
    return template


if __name__=="__main__":main()
