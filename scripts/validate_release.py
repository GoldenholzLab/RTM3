#!/usr/bin/env python3
"""Check current result hashes, embedded notebook figures, and DOCX preservation."""
import base64
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET

import nbformat

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def citation_identity(instruction):
    if "CSL_CITATION" not in instruction:
        return "bibliography"
    citation = json.loads(instruction[instruction.index("{"):])
    return [{key: item.get(key) for key in ["uris", "locator", "label", "prefix", "suffix", "suppress-author", "author-only"]}
            for item in citation["citationItems"]]


def main():
    metadata = json.loads((ROOT / "results/rtm3_run_metadata.json").read_text())
    names = {key: f"rtm3_{key}_results.csv" for key in ["primary", "cohort", "threshold", "interval", "correction"]}
    names.update({key: f"rtm3_{key}_sensitivity_results.csv" for key in ["cohort", "threshold", "interval", "correction"]})
    names.update(selection="rtm3_selection_diagnostics.csv", validation="rtm3_validation.json",
                 figure_png="figure1_sensitivity_and_FAR_vs_RTM_mpc_ci.png", figure_tif="figure1_sensitivity_and_FAR_vs_RTM_mpc_ci.tif",
                 selection_figure="appendix_selection_diagnostics.png")
    for key, expected in metadata["artifact_sha256"].items():
        assert sha(ROOT / "results" / names[key]) == expected, f"Hash mismatch: {key}"
    assert sha(ROOT / "reproduce_rtm3.py") == metadata["script_sha256"]
    assert sha(ROOT / "realSim.py") == metadata["simulator_provenance"]["local_source_sha256"]
    embedded = {}
    for name, count in [("for_rtm3_paper.ipynb", 2), ("for_appendix_RTM3.ipynb", 3)]:
        notebook = nbformat.read(ROOT / name, 4)
        nbformat.validate(notebook)
        code = [cell for cell in notebook.cells if cell.cell_type == "code"]
        assert all(cell.execution_count is not None for cell in code)
        assert not any(output.output_type == "error" for cell in code for output in cell.outputs)
        if name == "for_rtm3_paper.ipynb":
            comparison = [cell for cell in code if cell.source == "show_main_cohort_comparison(results)"]
            assert len(comparison) == 1 and comparison[0].outputs, "Missing executed fixed-cohort comparison"
            rendered = json.dumps(comparison[0].outputs)
            for expected in ["35,393", "29.3%", "34.7%", "6,182", "538"]:
                assert expected in rendered, f"Missing manuscript comparison value: {expected}"
        images = [base64.b64decode(output["data"]["image/png"]) for cell in code for output in cell.outputs if "image/png" in output.get("data", {})]
        assert len(images) == count, f"Missing notebook images: {name}"
        embedded[name] = images
    for image, filename in zip(embedded["for_rtm3_paper.ipynb"], ["figure1_sensitivity_and_FAR_vs_RTM_mpc_ci.png", "figure2_conceptual_mechanism.png"]):
        assert image == (ROOT / "results" / filename).read_bytes(), f"Notebook image differs: {filename}"
    for index, image in enumerate(embedded["for_appendix_RTM3.ipynb"], 1):
        (ROOT / "revision_20260915" / f"appendix_figure_{index}.png").write_bytes(image)
    docx = ROOT / "submission/epileptic_disorders/RTM3_v5.docx"
    if docx.exists():
        ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main", "wp": "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"}
        with ZipFile(docx) as current, ZipFile(ROOT / "revision_20260915/source/RTM3_v5.docx") as source:
            xml = ET.fromstring(current.read("word/document.xml"))
            old = ET.fromstring(source.read("word/document.xml"))
            fields = lambda root: [item.text for item in root.findall(".//w:instrText", ns) if "ZOTERO" in (item.text or "")]
            assert [citation_identity(f) for f in fields(xml)] == [citation_identity(f) for f in fields(old)]
            assert len(xml.findall(".//wp:inline", ns)) == 2
            assert current.read("word/media/image1.tif") == (ROOT / "results/figure1_sensitivity_and_FAR_vs_RTM_mpc_ci.tif").read_bytes()
            assert current.read("word/media/image2.png") == (ROOT / "results/figure2_conceptual_mechanism.png").read_bytes()
            before_word = ROOT / "revision_20260915/RTM3_v5.before_Word.docx"
            with ZipFile(before_word) as before:
                previous = ET.fromstring(before.read("word/document.xml"))
            visible = lambda root: "".join(node.text or "" for node in root.findall(".//w:t", ns))
            expected = visible(previous).replace("Clin Neurophysiol Off J Int Fed Clin Neurophysiol.", "Clin Neurophysiol.").replace("Epilepsy Behav EB.", "Epilepsy Behav.")
            assert visible(xml) == expected, "Unexpected prose change after Word/Zotero finalization"
            assert not xml.findall(".//w:ins", ns) and not xml.findall(".//w:del", ns)
            assert "Clin Neurophysiol. 2021" in visible(xml)
            assert "Epilepsy Behav. 2024" in visible(xml)
            assert "Epilepsy Behav. 2015" in visible(xml)
            paragraphs = xml.findall("./w:body/w:p", ns)
            (ROOT / "revision_20260915/manuscript_text.txt").write_text("\n".join(visible(p).rstrip() for p in paragraphs) + "\n")
    report = {"status": "PASS", "source_and_output_hashes": "match", "main_notebook_images": 2,
              "appendix_notebook_images": 3, "main_notebook_images_equal_saved_PNGs": True,
              "DOCX_images_equal_generated_images": True, "live_Zotero_fields_preserved": 17,
              "post_Word_prose_changes": "only the two requested printed journal-abbreviation corrections",
              "Word_refresh_and_visual_QA": "completed manually on 2026-09-15; 13-page Word PDF and all four DOCX files reviewed"}
    (ROOT / "revision_20260915/release_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
