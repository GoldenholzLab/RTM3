#!/usr/bin/env python3
import json
import random
import shutil
import string
import tempfile
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

DOCX = Path("epilepsia_brief_communication_draft.docx")
USER_ID = "1538114"

NS = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
}
ET.register_namespace("w", NS["w"])

W = f"{{{NS['w']}}}"


ITEMS = {
    1: (638, "2E3DNX5V", "Automated seizure detection using wearable devices: A clinical practice guideline of the International League Against Epilepsy and the International Federation of Clinical Neurophysiology"),
    2: (640, "WGUAU7CR", "Automated seizure detection with noninvasive wearable devices: A systematic review and meta-analysis"),
    3: (642, "DIR56C2D", "Quantifying and controlling the impact of regression to the mean on randomized controlled trials in epilepsy"),
    4: (644, "AGKVQQIK", "Response to placebo in clinical epilepsy trials--Old ideas and new insights"),
    5: (646, "9BM8NJ4B", "Confusing placebo effect with natural history in epilepsy: A big data approach"),
    6: (648, "EEUJRT6Z", "A multi-dataset time-reversal approach to clinical trial placebo response and the relationship to natural variability in epilepsy"),
    7: (650, "ZRBV3KWP", "Natural variability in seizure frequency: Implications for trials and placebo"),
    8: (652, "4S99MZFF", "Flexible realistic simulation of seizure occurrence recapitulating statistical properties of seizure diaries"),
    9: (654, "BK98X8PE", "Randomized phase 2 study of adjunctive cenobamate in patients with uncontrolled focal seizures"),
    10: (656, "VCKFKX24", "Safety and efficacy of adjunctive cenobamate (YKP3089) in patients with uncontrolled focal seizures: a multicentre, double-blind, randomised, placebo-controlled, dose-response trial"),
    11: (658, "QQJQYVX2", "Simulating Clinical Trials With and Without Intracranial EEG Data"),
    12: (661, "AA8BTYPI", "Factors determining response to antiepileptic drugs in randomized controlled trials. A systematic review and meta-analysis"),
    13: (663, "EJCH3TFJ", "Placebo and nocebo responses in drug trials of epilepsy"),
    14: (665, "8DQ4Y9T9", "How much of the placebo 'effect' is really statistical regression?"),
}


CITATION_SEQUENCE = [
    ("1,2", [1, 2]),
    ("3-7", [3, 4, 5, 6, 7]),
    ("3", [3]),
    ("5-7", [5, 6, 7]),
    ("8", [8]),
    ("9,10", [9, 10]),
    ("3", [3]),
    ("11", [11]),
    ("1,2", [1, 2]),
    ("8", [8]),
]


def rid(prefix="z"):
    return prefix + "".join(random.choice(string.ascii_letters + string.digits) for _ in range(8))


def item_payload(ref_num):
    item_id, key, title = ITEMS[ref_num]
    return {
        "id": item_id,
        "uris": [f"http://zotero.org/users/{USER_ID}/items/{key}"],
        "itemData": {
            "id": item_id,
            "type": "article-journal",
            "title": title,
        },
    }


def citation_json(ref_nums, formatted):
    return json.dumps(
        {
            "citationID": rid("c"),
            "properties": {
                "formattedCitation": f"<sup>{formatted}</sup>",
                "plainCitation": formatted,
                "noteIndex": 0,
            },
            "citationItems": [item_payload(n) for n in ref_nums],
            "schema": "https://github.com/citation-style-language/schema/raw/master/csl-citation.json",
        },
        ensure_ascii=False,
        separators=(",", ":"),
    )


def bibliography_json():
    cited = {n for _, nums in CITATION_SEQUENCE for n in nums}
    uncited = [
        [f"http://zotero.org/users/{USER_ID}/items/{ITEMS[n][1]}"]
        for n in sorted(set(ITEMS) - cited)
    ]
    return json.dumps(
        {"uncited": uncited, "omitted": [], "custom": []},
        ensure_ascii=False,
        separators=(",", ":"),
    )


def paragraph_text(p):
    return "".join(t.text or "" for t in p.findall(".//w:t", NS))


def is_superscript_run(r):
    va = r.find("w:rPr/w:vertAlign", NS)
    return va is not None and va.get(W + "val") == "superscript"


def make_text_run(text, superscript=False):
    r = ET.Element(W + "r")
    if superscript:
        rpr = ET.SubElement(r, W + "rPr")
        ET.SubElement(rpr, W + "vertAlign", {W + "val": "superscript"})
    t = ET.SubElement(r, W + "t")
    if text.startswith(" ") or text.endswith(" "):
        t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    t.text = text
    return r


def make_field_runs(zotero_code, visible_text="", superscript=False):
    begin = ET.Element(W + "r")
    ET.SubElement(begin, W + "fldChar", {W + "fldCharType": "begin"})

    instr = ET.Element(W + "r")
    instr_text = ET.SubElement(instr, W + "instrText")
    instr_text.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    instr_text.text = f" ADDIN ZOTERO_{zotero_code} "

    separate = ET.Element(W + "r")
    ET.SubElement(separate, W + "fldChar", {W + "fldCharType": "separate"})

    result = make_text_run(visible_text, superscript=superscript)

    end = ET.Element(W + "r")
    ET.SubElement(end, W + "fldChar", {W + "fldCharType": "end"})

    return [begin, instr, separate, result, end]


def replace_run_with_field(p, run, field_runs):
    children = list(p)
    idx = children.index(run)
    p.remove(run)
    for offset, new_run in enumerate(field_runs):
        p.insert(idx + offset, new_run)


def insert_doc_var(settings_root):
    doc_vars = settings_root.find("w:docVars", NS)
    if doc_vars is None:
        doc_vars = ET.SubElement(settings_root, W + "docVars")
    for existing in list(doc_vars):
        if existing.get(W + "name") == "ZOTERO_PREF":
            doc_vars.remove(existing)
    value = (
        '<data data-version="3" zotero-version="7.0">'
        f'<session id="{rid("s")}"/>'
        '<style id="http://www.zotero.org/styles/epilepsia" hasBibliography="1" bibliographyStyleHasBeenSet="1"/>'
        '<prefs>'
        '<pref name="fieldType" value="Field"/>'
        '<pref name="noteType" value="0"/>'
        '<pref name="automaticJournalAbbreviations" value="true"/>'
        "</prefs></data>"
    )
    ET.SubElement(doc_vars, W + "docVar", {W + "name": "ZOTERO_PREF", W + "val": value})


def main():
    if not DOCX.exists():
        raise SystemExit(f"Missing {DOCX}")

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        with zipfile.ZipFile(DOCX) as zin:
            zin.extractall(tmpdir)

        doc_path = tmpdir / "word" / "document.xml"
        settings_path = tmpdir / "word" / "settings.xml"

        doc_tree = ET.parse(doc_path)
        root = doc_tree.getroot()
        body = root.find("w:body", NS)
        paragraphs = body.findall("w:p", NS)

        intro_idx = next(i for i, p in enumerate(paragraphs) if paragraph_text(p).strip() == "Introduction")
        refs_idx = next(i for i, p in enumerate(paragraphs) if paragraph_text(p).strip() == "References")

        search_start = intro_idx + 1
        for visible, refs in CITATION_SEQUENCE:
            found = False
            for p in paragraphs[search_start:refs_idx]:
                for r in list(p.findall("w:r", NS)):
                    if not is_superscript_run(r):
                        continue
                    if paragraph_text(r) != visible:
                        continue
                    code = "ITEM CSL_CITATION " + citation_json(refs, visible)
                    replace_run_with_field(p, r, make_field_runs(code, visible, superscript=True))
                    search_start = paragraphs.index(p)
                    found = True
                    break
                if found:
                    break
            if not found:
                raise RuntimeError(f"Could not find superscript citation {visible!r}")

        fig_idx = next(i for i, p in enumerate(paragraphs[refs_idx + 1 :], refs_idx + 1) if "Figure 1." in paragraph_text(p))
        for p in paragraphs[refs_idx + 1 : fig_idx]:
            body.remove(p)

        bib_p = ET.Element(W + "p")
        bib_code = "BIBL " + bibliography_json() + " CSL_BIBLIOGRAPHY"
        placeholder = "Zotero bibliography will be generated on refresh."
        for r in make_field_runs(bib_code, placeholder, superscript=False):
            bib_p.append(r)
        body.insert(list(body).index(paragraphs[refs_idx]) + 1, bib_p)

        settings_tree = ET.parse(settings_path)
        insert_doc_var(settings_tree.getroot())

        doc_tree.write(doc_path, encoding="UTF-8", xml_declaration=True, short_empty_elements=True)
        settings_tree.write(settings_path, encoding="UTF-8", xml_declaration=True, short_empty_elements=True)

        new_docx = tmpdir / DOCX.name
        with zipfile.ZipFile(new_docx, "w", compression=zipfile.ZIP_DEFLATED) as zout:
            for path in tmpdir.rglob("*"):
                if path == new_docx or path.is_dir():
                    continue
                zout.write(path, path.relative_to(tmpdir))
        shutil.move(new_docx, DOCX)


if __name__ == "__main__":
    main()
