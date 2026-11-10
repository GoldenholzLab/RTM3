#!/usr/bin/env python3
"""Validate presentation outputs against the pinned release without simulation."""
import ast
import base64
import csv
import hashlib
import json
import math
import subprocess
import types
from fractions import Fraction
from decimal import Decimal
from pathlib import Path
from zipfile import ZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches
from lxml import etree

from export_appendix import BASE_COMMIT, NS, ROOT, table_rows

OUT=ROOT/"results"/"appendix_export"
WARNING="The fixed perfect-detector cohort is hypothetical and completely unrealistic as an enrollment strategy for a trial whose eligibility is determined by the imperfect detector being evaluated."


def blob(path):
    return subprocess.check_output(["git","show",f"{BASE_COMMIT}:{path}"],cwd=ROOT)


def sha(data):return hashlib.sha256(data).hexdigest()
def norm(value):return " ".join(value.replace("\u200b", "").split())


def verify_tables(nb, original):
    def tables(book):
        return [table_rows(o.data["text/html"]) for c in book.cells for o in c.get("outputs",[]) if "text/html" in o.get("data",{})]
    old=tables(original);new=tables(nb)
    assert len(old)==21 and len(new)==22
    unmatched=list(range(len(new)))
    coverage=[]
    for i,table in enumerate(old,1):
        match=next(j for j in unmatched if new[j]==table)
        unmatched.remove(match)
        coverage.append({"original_table":i,"updated_table":match+1,"body_rows":len(table[1]),"all_cells_identical":True})
    assert len(unmatched)==1 and "Retained n" in new[unmatched[0]][0]
    doc=Document(OUT/"RTM3_appendix.docx")
    assert len(doc.tables)==22
    for i,(table,expected) in enumerate(zip(doc.tables,new),1):
        actual=[[norm(c.text) for c in r.cells] for r in table.rows]
        want=[[norm(x) for x in expected[0]]]+[[norm(x) for x in r] for r in expected[1]]
        assert actual==want,(i,actual[:2],want[:2])
        assert table.rows[0]._tr.xpath("./w:trPr/w:tblHeader")
        for row_index, row in enumerate(table.rows):
            assert row._tr.xpath("./w:trPr/w:cantSplit")
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    assert paragraph.paragraph_format.keep_with_next == (row_index == 0)
    html=BeautifulSoup((ROOT/"results/notebooks/for_appendix_RTM3.html").read_text(),"html.parser")
    assert len(html.select("table"))==22
    for tag,expected in zip(html.select("table"),new):
        assert table_rows(str(tag))==expected
    return {"original_tables":21,"retained_original_body_rows":sum(len(t[1]) for t in old),
            "updated_tables":22,"updated_body_rows":sum(len(t[1]) for t in new),"coverage":coverage}


def verify_figures(nb):
    import sys
    sys.path.insert(0,str(ROOT))
    import rtm3_reporting as revised
    original=types.ModuleType("original_reporting")
    exec(compile(blob("rtm3_reporting.py"),"original_reporting.py","exec"),original.__dict__)
    results=revised.load_results(ROOT/"results")
    configurations=[(results["cohort"],"cohort_strategy","Fixed versus reselected cohorts"),
                    (results["correction"],"correction_strategy","False-alarm correction and interval alternatives")]
    orig_rule=results["interval"].query("eligibility_variant == 'Original 25-day rule retained'")
    configurations.append((pd.concat([results["primary"].assign(rule="Main: monthly counts only"),orig_rule.assign(rule="Monthly counts plus original 25-day rule")]),"rule","Sensitivity analysis of the 25-day eligibility rule"))
    for frame,group,title in configurations:
        a=original.comparison_figure(frame,group,title);b=revised.comparison_figure(frame,group,title)
        for i,(aa,bb) in enumerate(zip(a.axes,b.axes)):
            assert len(aa.lines)==len(bb.lines)
            for la,lb in zip(aa.lines,bb.lines):
                np.testing.assert_array_equal(la.get_xdata(),lb.get_xdata());np.testing.assert_array_equal(la.get_ydata(),lb.get_ydata())
                assert (la.get_label(),la.get_marker(),la.get_linestyle())==(lb.get_label(),lb.get_marker(),lb.get_linestyle())
            assert aa.get_xlim()==bb.get_xlim()
            if group=="cohort_strategy" and i==3:
                assert bb.get_ylim()==(0,40)
                np.testing.assert_array_equal(bb.get_yticks(),[0,10,20,30,40])
            else:
                assert aa.get_ylim()==bb.get_ylim()
                np.testing.assert_array_equal(aa.get_yticks(),bb.get_yticks())
        plt.close(a);plt.close(b)
    pngs=[base64.b64decode(o.data["image/png"]) for c in nb.cells for o in c.get("outputs",[]) if "image/png" in o.get("data",{})]
    assert len(pngs)==3
    html=BeautifulSoup((ROOT/"results/notebooks/for_appendix_RTM3.html").read_text(),"html.parser")
    html_pngs=[base64.b64decode(i["src"].split(",",1)[1]) for i in html.select("img")]
    assert [sha(i) for i in pngs]==[sha(i) for i in html_pngs]
    with ZipFile(OUT/"RTM3_appendix.docx") as z:
        word_hashes={sha(z.read(n)) for n in z.namelist() if n.startswith("word/media/")}
    assert word_hashes=={sha(i) for i in pngs}
    doc=Document(OUT/"RTM3_appendix.docx")
    assert len(doc.inline_shapes)==3
    assert all(shape.width==Inches(9.4) for shape in doc.inline_shapes)
    assert all(p.alignment==WD_ALIGN_PARAGRAPH.CENTER for p in doc.paragraphs if p._p.xpath(".//w:drawing"))
    assert pngs[0]==(ROOT/"results/appendix_fixed_vs_reselected.png").read_bytes()
    assert not html.select("script[src],link[href]")
    assert all(not i["src"].startswith(("http:","https:","file:")) for i in html.select("img"))
    assert len(html.select("math"))==51
    return {"count":3,"all_plot_values_unchanged":True,"only_axis_change":"cohort FAR selection: 0 to 40; ticks 0,10,20,30,40",
            "PNG_bytes_identical_in_notebook_HTML_DOCX":True,"corrected_figure_sha256":sha(pngs[0]),"HTML_external_resources":0,
            "Word_inline_figure_width_inches":9.4,"Word_drawing_paragraphs_centered":True}


def verify_numbers():
    files=list((ROOT/"results").glob("*.csv"))+list((ROOT/"results").glob("*.json"))
    protected={str(p.relative_to(ROOT)):sha(p.read_bytes()) for p in files}
    for name in ["reproduce_rtm3.py","realSim.py"]:
        protected[name]=sha((ROOT/name).read_bytes())
    for path,digest in protected.items():assert digest==sha(blob(path)),path
    primary=pd.read_csv(ROOT/"results/rtm3_primary_results.csv")
    for p in (ROOT/"results").glob("*.csv"):
        frame=pd.read_csv(p)
        assert (frame.n_valid_MPC+frame.n_zero_baseline).eq(frame.n_eligible).all()
        assert frame.n_undefined_MPC.eq(frame.n_zero_baseline).all()
        np.testing.assert_allclose(frame.eligibility_fraction,frame.n_eligible/frame.n_total,rtol=0,atol=1e-15)
        if "n_observed_RTM" in frame:
            np.testing.assert_allclose(frame.observed_RTM_fraction,frame.n_observed_RTM/frame.n_eligible,rtol=0,atol=1e-15)
    cohort=pd.read_csv(ROOT/"results/rtm3_cohort_sensitivity_results.csv")
    fixed=cohort.query("cohort_strategy == 'Fixed perfect-detector cohort'")
    assert fixed.n_eligible.eq(35393).all()
    endpoints=[]
    for strategy in cohort.cohort_strategy.unique():
        part=cohort[cohort.cohort_strategy.eq(strategy)]
        rows=[part.query("sweep=='FAR' and FAR==0").iloc[0],part.query("sweep=='sensitivity' and sensitivity==.1").iloc[0],part.query("sweep=='FAR' and FAR==1").iloc[0]]
        expected=[(37.2,14.3,0),(29.3,33.3,6182),(34.7,10.4,538)] if strategy.startswith("Fixed") else [(37.2,14.3,0),(64.4,48.1,0),(53.6,31.2,0)]
        for row,want in zip(rows,expected):
            got=(round(100*row.observed_RTM_fraction,1),round(row.MPC_median,1),int(row.n_zero_baseline));assert got==want
            endpoints.append({"strategy":strategy,"sensitivity":row.sensitivity,"FAR":row.FAR,"observed_RTM_percent":got[0],"MPC_percent":got[1],"selected":int(row.n_eligible),"MPC_denominator":int(row.n_valid_MPC),"zero_baseline":got[2],"retained_failing_eligibility":int(row.n_selected_failing_current_eligibility)})
    thresholds=pd.read_csv(ROOT/"results/rtm3_threshold_sensitivity_results.csv")
    assert thresholds.n_eligibility_mask_disagreements.eq(0).all()
    a,b=[g.sort_values("FAR") for _,g in thresholds.groupby("threshold_strategy")]
    common=["n_eligible","n_observed_RTM","observed_RTM_fraction","MPC_median","MPC_ci_lo","MPC_ci_hi","n_valid_MPC"]
    np.testing.assert_array_equal(a[common].to_numpy(),b[common].to_numpy())
    corrections=pd.read_csv(ROOT/"results/rtm3_correction_sensitivity_results.csv")
    oracle=corrections[corrections.correction_strategy.eq("Exact oracle false-alarm removal")]
    assert all(oracle[c].nunique()==1 for c in common)
    original=pd.read_csv(ROOT/"results/rtm3_interval_sensitivity_results.csv").query("eligibility_variant=='Original 25-day rule retained'")
    paired=primary.merge(original,on=["sweep","sensitivity","FAR"],validate="one_to_one",suffixes=("_a","_b"))
    metadata=json.loads((ROOT/"results/rtm3_run_metadata.json").read_text())
    for c,m,scale in [("eligibility_fraction","eligibility",100),("observed_RTM_fraction","observed_RTM",100),("MPC_median","MPC",1)]:
        delta=(paired[c+"_a"]-paired[c+"_b"]).abs().max()*scale
        assert np.isclose(delta,metadata["validation"][f"max_absolute_{m}_difference_primary_vs_original_25day_pp"],rtol=0,atol=1e-12)
    selection=pd.read_csv(ROOT/"results/rtm3_selection_diagnostics.csv")
    for c in selection.columns:pd.testing.assert_series_equal(selection[c],primary[c],check_names=False)
    # Execute only the released pure endpoint function on small boundary fixtures.
    tree=ast.parse((ROOT/"reproduce_rtm3.py").read_text())
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="observed_rtm_indicator")
    scope={"np":np,"math":math,"Fraction":Fraction,"LONG_TERM_MONTHS":36,"BASELINE_MONTHS":2,"TEST_MONTHS":3}
    exec(compile(ast.Module(body=[fn],type_ignores=[]),"released_endpoint_only","exec"),scope)
    f=scope["observed_rtm_indicator"]
    monthly=np.array([[6,6,4,4,4],[6,6,5,5,5],[5,5,5,5,5],[6,6,6,6,6],[6,6,7,7,7]])
    assert f(monthly,np.full(5,5.),1.,long_term_total=np.full(5,180)).tolist()==[False,True,False,False,False]
    assert f(np.array([[1,1,0,0,0],[1,1,1,0,0]]),np.full(2,5.),.1,long_term_total=np.full(2,180)).tolist()==[False,True]
    return {"protected_file_hashes":protected,"protected_files_byte_identical_to_pinned_release":True,"endpoint_values":endpoints,"CSV_denominators_and_ratios":True,"threshold_equivalence":True,"oracle_invariance":True,"interval_maxima_reconciled":True,"selection_CSV_agreement":True,"strict_tie_boundary_fixtures":7,"simulation_rerun":False}


def verify_interval_grid_comparison():
    """Independently check the comparison grids using CSV counts and decimals."""
    with (ROOT/"results/rtm3_primary_results.csv").open() as f:
        primary=list(csv.DictReader(f))
    with (ROOT/"results/rtm3_interval_sensitivity_results.csv").open() as f:
        interval=list(csv.DictReader(f))
    def key(row):return (row["sweep"],Decimal(row["sensitivity"]),Decimal(row["FAR"]))
    lookup={key(row):row for row in primary}
    def maxima(rows):
        rtm=[];mpc=[]
        for row in rows:
            reference=lookup[key(row)]
            # Derive prevalence from integer numerator/denominator, independently
            # of the stored fractions used by the reporting tables.
            r=abs(Fraction(int(row["n_observed_RTM"]),int(row["n_eligible"]))-
                  Fraction(int(reference["n_observed_RTM"]),int(reference["n_eligible"])))*100
            m=abs(Decimal(row["MPC_median"])-Decimal(reference["MPC_median"]))
            rtm.append((r,key(row)));mpc.append((m,key(row)))
        answer={"n_conditions":len(rows)}
        for metric,values in [("observed_RTM",rtm),("MPC",mpc)]:
            maximum,condition=max(values,key=lambda pair:pair[0])
            answer[metric]={"maximum_absolute_change_pp":float(maximum),
                            "formatted_10dp":f"{float(maximum):.10f}",
                            "sweep":condition[0],"sensitivity":str(condition[1]),"FAR":str(condition[2])}
        return answer
    original=[r for r in interval if r["eligibility_variant"]=="Original 25-day rule retained"]
    observed=[r for r in interval if r["eligibility_variant"]=="25-day rule on observed days"]
    oracle=[r for r in interval if r["eligibility_variant"]=="25-day rule on latent true days"]
    assert len(original)==21 and len(observed)==len(oracle)==11
    assert all(r["sweep"]=="FAR" for r in observed+oracle)
    common=maxima([r for r in original if r["sweep"]=="FAR"])
    observed_max=maxima(observed);oracle_max=maxima(oracle);full=maxima(original)
    for metric,expected in [("observed_RTM","0.9571218755"),("MPC","0.8658008658")]:
        assert common[metric]["formatted_10dp"]==observed_max[metric]["formatted_10dp"]==expected
        assert common[metric]["maximum_absolute_change_pp"]==observed_max[metric]["maximum_absolute_change_pp"]
        assert full[metric]["sweep"]=="sensitivity"
        assert oracle_max[metric]["maximum_absolute_change_pp"]>common[metric]["maximum_absolute_change_pp"]
    return {"method":"RTM prevalence recalculated from integer counts with exact Fraction arithmetic; MPC differences from Decimal CSV strings",
            "original_rule_full_grid":full,"original_rule_common_FAR_grid":common,
            "observed_day_common_FAR_grid":observed_max,"oracle_day_FAR_grid":oracle_max,
            "common_grid_maxima_equal":True,"alternatives_only_evaluated_on_11_point_FAR_grid":True}


def verify_structure(nb,original):
    with ZipFile(OUT/"RTM3_appendix.docx") as z:
        root=etree.fromstring(z.read("word/document.xml"))
        for name in z.namelist():
            if name.endswith(".xml"):etree.fromstring(z.read(name))
        text="".join(root.xpath("//w:t/text()",namespaces=NS))
        assert WARNING in text
        wp=root.xpath("//w:p[.//w:t[contains(.,'The fixed perfect-detector cohort is hypothetical')]]",namespaces=NS)
        assert len(wp)==2
        for p in wp:
            for r in p.xpath("./w:r[w:t]",namespaces=NS):assert r.xpath("./w:rPr/w:b",namespaces=NS)
        math=root.xpath("//m:oMath",namespaces=NS);assert len(math)==51
        seq=root.xpath("//w:instrText[contains(.,'SEQ Equation')]",namespaces=NS);assert len(seq)==9
        cache=[]
        for q in seq:
            p=q.getparent().getparent();cache.append("".join(p.xpath("./w:r/w:t/text()",namespaces=NS)).strip())
        assert cache==[f"({i})" for i in range(1,10)],cache
        eqs=json.loads((OUT/"native_equations/equations.json").read_text())
        displays={e["live_equation_number"]:e for e in eqs if e["kind"]=="DisplayMath"}
        assert displays[3]["omml_text"].count("≥")==3
        assert displays[4]["omml_text"].count("≥")==3
        assert displays[6]["omml_text"].count(">")==2
        endpoint=etree.parse(str(OUT/"native_equations"/displays[6]["omml_file"]))
        assert endpoint.xpath("//m:dPr/m:begChr[@m:val='|']",namespaces=NS)
        assert endpoint.xpath("//m:dPr/m:endChr[@m:val='|']",namespaces=NS)
        assert len(root.xpath("//w:tbl/w:tr[1]/w:trPr/w:tblHeader",namespaces=NS))==22
        # No prior DOCX or Zotero fields were imported; no fields were flattened.
        allfields=root.xpath("//w:instrText/text()",namespaces=NS)
        assert not any("ZOTERO" in x for x in allfields)
    old_md=[o.data["text/markdown"] for c in original.cells for o in c.get("outputs",[]) if "text/markdown" in o.get("data",{})]
    new_md=[o.data["text/markdown"] for c in nb.cells for o in c.get("outputs",[]) if "text/markdown" in o.get("data",{})]
    old_grid_comparison = 'Compared with the main analysis, the observed-day rule changed observed RTM and MPC less than the original rule did; the oracle-day rule changed them more. These comparisons concern the largest absolute differences across the sweep, not a uniform ordering at every condition.'
    revised_grid_comparison = "The original rule's larger overall maxima include the sensitivity sweep; the observed-day and oracle-day alternatives were evaluated only on the 11-point FAR sweep. On that common FAR sweep, the observed-day and original rules have equal maximum absolute changes from the main analysis: 0.9571218755 percentage points for observed RTM and 0.8658008658 percentage points for MPC. The oracle-day rule has larger maxima on that same FAR grid. These comparisons concern maxima across conditions, not a uniform ordering at every condition."
    assert sum(old_grid_comparison in x for x in old_md) == 1
    assert all(x.replace(old_grid_comparison, revised_grid_comparison) in new_md for x in old_md)
    old_plain=[o.data["text/plain"] for c in original.cells for o in c.get("outputs",[]) if set(o.get("data",{}))=={"text/plain"}]
    new_plain=[o.data["text/plain"] for c in nb.cells for o in c.get("outputs",[]) if set(o.get("data",{}))=={"text/plain"}]
    assert old_plain==new_plain
    assert sum(WARNING in c.source for c in nb.cells if c.cell_type=="markdown")==2
    html=BeautifulSoup((ROOT/"results/notebooks/for_appendix_RTM3.html").read_text(),"html.parser")
    assert sum(norm(t.get_text())==WARNING for t in html.select("strong"))==2
    return {"OOXML_well_formed":True,"bold_warning_in_three_formats":True,"original_generated_narratives_accounted_for":len(old_md),"intentional_grid_qualification_corrections":1,"original_plain_provenance_validation_outputs_retained":len(old_plain),"equations_native":51,"display_equations_with_live_numbers":9,"strict_and_inclusive_inequalities_and_absolute_value_structure":True,"repeated_table_headers":22,"Zotero_fields_in_source":"None; inputs were Markdown, HTML, notebook, CSV, JSON"}


def main():
    nb=nbformat.read(ROOT/"for_appendix_RTM3.ipynb",as_version=4);nbformat.validate(nb)
    original=nbformat.reads(blob("for_appendix_RTM3.ipynb").decode(),as_version=4)
    audit={"status":"PASS structural and content checks; pending authoritative Word verification",
           "tables":verify_tables(nb,original),"figures":verify_figures(nb),"numbers":verify_numbers(),"structure":verify_structure(nb,original),
           "interval_grid_comparison":verify_interval_grid_comparison()}
    (OUT/"validation_audit.json").write_text(json.dumps(audit,indent=2))
    print(json.dumps({k:v if k=="status" else "PASS" for k,v in audit.items()},indent=2))


if __name__=="__main__":main()
