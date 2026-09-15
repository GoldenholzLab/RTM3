#!/usr/bin/env python3
"""Apply layout-only repairs to closed Word submission documents.

Run after Zotero refresh. Reopen, save, and export through Word afterward.
The script preserves all text, citation instructions, media, and relationships.
"""
import copy
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

from lxml import etree

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "submission/epileptic_disorders/RTM3_v5.docx"
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def text(paragraph):
    return "".join(x.text or "" for x in paragraph.iter(W + "t"))


def properties(parent, name):
    element = parent.find(W + name)
    if element is None:
        element = etree.Element(W + name)
        parent.insert(0, element)
    return element


def setting(parent, name, **attributes):
    element = parent.find(W + name)
    if element is None:
        element = etree.SubElement(parent, W + name)
    for key, value in attributes.items():
        element.set(W + key, str(value))
    return element


def label_only_bold(paragraph, label):
    # These paragraphs contain plain text and, for Figure 1, a preserved drawing.
    for run in list(paragraph.findall(W + "r")):
        element = run.find(W + "t")
        if element is None:
            continue
        value = element.text or ""
        setting(properties(run, "rPr"), "b", val="0")
        if value.startswith(label):
            label_run = copy.deepcopy(run)
            label_run.find(W + "t").text = label
            setting(properties(label_run, "rPr"), "b", val="1")
            element.text = value[len(label):]
            element.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
            run.addprevious(label_run)


def main():
    with ZipFile(PATH) as archive:
        entries = [(entry, archive.read(entry.filename)) for entry in archive.infolist()]
    root = etree.fromstring(dict((entry.filename, data) for entry, data in entries)["word/document.xml"])
    original_text = text(root)
    original_fields = [node.text for node in root.iter(W + "instrText")]
    paragraphs = root.find(W + "body").findall(W + "p")
    in_education = False
    in_key_points = False
    for paragraph in paragraphs:
        value = text(paragraph)
        for label in ["Objective:", "Methods:", "Results:", "Significance:", "Figure 1.", "Figure 2."]:
            if value.startswith(label):
                label_only_bold(paragraph, label)
        if value == "Key points":
            in_key_points = True
        elif in_key_points and not value.startswith("• "):
            in_key_points = False
        if in_key_points:
            ppr = properties(paragraph, "pPr")
            border = setting(ppr, "pBdr")
            for edge in ["top", "left", "bottom", "right"]:
                setting(border, edge, val="single", sz="4", space="6", color="808080")
            setting(ppr, "shd", val="clear", fill="F2F2F2")
        if value == "Test Yourself":
            in_education = True
            setting(properties(paragraph, "pPr"), "pageBreakBefore", val="1")
        if in_education:
            ppr = properties(paragraph, "pPr")
            setting(ppr, "spacing", line="240", lineRule="auto", after="120")
            setting(ppr, "keepLines", val="1")
            setting(ppr, "keepNext", val="0" if value.startswith("Answer:") else "1")
    assert text(root) == original_text
    assert [node.text for node in root.iter(W + "instrText")] == original_fields
    document = etree.tostring(root, encoding="UTF-8", xml_declaration=True, standalone=True)
    with ZipFile(PATH, "w", ZIP_DEFLATED) as archive:
        for entry, data in entries:
            archive.writestr(entry, document if entry.filename == "word/document.xml" else data)
    print("Layout repaired; prose, fields, media, and relationships preserved.")


def polish_title_page():
    path = PATH.with_name("titlepage.docx")
    with ZipFile(path) as archive:
        entries = [(entry, archive.read(entry.filename)) for entry in archive.infolist()]
    root = etree.fromstring(dict((entry.filename, data) for entry, data in entries)["word/document.xml"])
    original_text = text(root)
    for index, paragraph in enumerate(root.find(W + "body").findall(W + "p")):
        ppr = properties(paragraph, "pPr")
        setting(ppr, "spacing", line="240", lineRule="auto", after="120")
        setting(ppr, "jc", val="left")
        for run in paragraph.findall(W + "r"):
            rpr = properties(run, "rPr")
            setting(rpr, "b", val="1" if index == 0 or text(paragraph) in ["Affiliations:", "Author emails"] else "0")
            setting(rpr, "sz", val="28" if index == 0 else "24")
    assert text(root) == original_text
    document = etree.tostring(root, encoding="UTF-8", xml_declaration=True, standalone=True)
    with ZipFile(path, "w", ZIP_DEFLATED) as archive:
        for entry, data in entries:
            archive.writestr(entry, document if entry.filename == "word/document.xml" else data)
    print("Title-page spacing repaired; author details and prose preserved.")


def polish_disclosure():
    path = PATH.with_name("Conflict of interest.docx")
    with ZipFile(PATH) as main_file:
        manuscript = etree.fromstring(main_file.read("word/document.xml"))
    statement = next(text(p) for p in manuscript.iter(W + "p") if text(p).startswith("DMG has been provided speaker fees"))
    statement = statement.replace("DMG has", "Daniel M. Goldenholz has").replace("grants from NIH, ABPN,", "grants from the National Institutes of Health, the American Board of Psychiatry and Neurology,")
    with ZipFile(path) as archive:
        entries = [(entry, archive.read(entry.filename)) for entry in archive.infolist()]
    root = etree.fromstring(dict((entry.filename, data) for entry, data in entries)["word/document.xml"])
    paragraph = next(p for p in root.iter(W + "p") if text(p).startswith("DMG has been provided speaker fees"))
    run = copy.deepcopy(paragraph.find(W + "r"))
    for child in list(run):
        if child.tag != W + "rPr":
            run.remove(child)
    etree.SubElement(run, W + "t").text = statement
    for child in list(paragraph):
        if child.tag != W + "pPr":
            paragraph.remove(child)
    paragraph.append(run)
    document = etree.tostring(root, encoding="UTF-8", xml_declaration=True, standalone=True)
    with ZipFile(path, "w", ZIP_DEFLATED) as archive:
        for entry, data in entries:
            archive.writestr(entry, document if entry.filename == "word/document.xml" else data)
    print("Standalone disclosure now matches the manuscript, with abbreviations expanded.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--title-only", action="store_true")
    parser.add_argument("--disclosure-only", action="store_true")
    args = parser.parse_args()
    if args.disclosure_only:
        polish_disclosure()
    elif args.title_only:
        polish_title_page()
    else:
        main()
        polish_title_page()
