#!/usr/bin/env python3
"""Put the bar plot I and J into the hand-edited figure, leaving panels A-H as edited.

Input is earlier_versions/Figure1_handedited_scatter_IJ.svg (Inkscape edits, scatter I and J),
so this can be rerun. The labels "11 h / 15 h after onset" had been moved under I and J; they go
back above G and H, because I and J now carry their own, and G and H move down slightly to make
room. Writes Figure1.svg and Figure1.png (2600 x 3100 px, the size of the hand export).

Also applies the label changes to the hand-edited text: Control Dawn/Dusk in A and B, sentence-case Dawn and Dusk in C-E, hours after cold onset in
C and D, DEGs in E. G and H are replaced with the current versions from the code layout, keeping
their hand-set position and size.
"""
import re
import subprocess

from lxml import etree

from figure1_style import DATA, ROOT

SVG = "{http://www.w3.org/2000/svg}"
SOURCE = ROOT / "earlier_versions" / "Figure1_handedited_scatter_IJ.svg"
BARS = DATA / "panel_I_J" / "Figure1I_J_bars.svg"
CODE_LAYOUT = ROOT / "earlier_versions" / "Figure1_code_layout_scatter_IJ.svg"
X_PT, Y_PT = 231, 400.5  # top-left of the inserted block, in the 468 x 558 pt page
G_H_SHIFT = 4.2  # pt down, clears "Time point (h)" of panel D
E_KEY_SHIFT = -4  # pt left, keeps the key of panel E clear of the moved axis of G
TITLES = (("11 h after onset", 298), ("15 h after onset", 414.4))  # centres of the G and H axes
TITLE_BASELINE = 272


def text_of(node):
    return "".join(node.itertext()).strip()


def prefixed(root, prefix):
    """Copy of root's content with every id renamed, so nothing clashes with the figure's ids."""
    raw = etree.tostring(root, encoding="unicode")
    raw = re.sub(r'\bid="([^"]+)"', rf'id="{prefix}\1"', raw)
    raw = re.sub(r'href="#([^"]+)"', rf'href="#{prefix}\1"', raw)
    raw = re.sub(r"url\(#([^)]+)\)", rf"url(#{prefix}\1)", raw)
    return etree.fromstring(raw)


def replace_text(node, old, new, expected):
    """Replace whole-string text anywhere under node; stop if the count is not as expected."""
    hits = [el for el in node.iter(f"{SVG}text", f"{SVG}tspan") if el.text and el.text.strip() == old]
    assert len(hits) == expected, (old, len(hits))
    for el in hits:
        el.text = new


def main(code_layout=CODE_LAYOUT, out=ROOT / "Figure1.svg"):
    tree = etree.parse(str(SOURCE))
    root = tree.getroot()
    figure = root.find(f"{SVG}g[@id='figure_1']")
    groups = {g.get("id"): g for g in figure.findall(f"{SVG}g")}

    for key in ("axes_11", "axes_12", "legend_6"):
        figure.remove(groups[key])
    letters = [g for k, g in groups.items() if k.startswith("text_") and text_of(g) in ("I", "J")]
    assert len(letters) == 2
    for g in letters:
        figure.remove(g)

    replace_text(groups["axes_1"], "Cntl Dawn", "Control Dawn", 1)
    replace_text(groups["axes_1"], "Cold Dawn d1", "Cold Dawn", 1)
    replace_text(groups["axes_3"], "Cntl Dusk", "Control Dusk", 1)
    replace_text(groups["axes_3"], "Cold Dusk d1", "Cold Dusk", 1)
    for axes, key in (("axes_1", "legend_1"), ("axes_3", "legend_2")):  # Control is wider than Cntl
        entries = groups[axes].find(f"{SVG}g[@id='{key}']").findall(f"{SVG}g")
        start = next(i for i, g in enumerate(entries) if "Control" in text_of(g)) + 1
        assert len(entries) - start == 2, key
        for g in entries[start:]:
            g.set("transform", f"translate(11,0) {g.get('transform', '')}".strip())
    for old, new in (("DAWN vs CONT", "Dawn vs Control"), ("DUSK vs CONT", "Dusk vs Control"), ("DAWN up", "Dawn up"),
                     ("DAWN down", "Dawn down"), ("DUSK up", "Dusk up"), ("DUSK down", "Dusk down")):
        replace_text(figure, old, new, 1)
    replace_text(figure, "Time point (h)", "Hours after cold onset", 2)
    replace_text(figure, "Number of DEGs", "DEGs", 1)
    d_axis = groups["axes_6"].find(f"{SVG}g[@id='matplotlib.axis_11']")
    for hour in range(15, 40, 4):
        replace_text(d_axis, str(hour), str(hour - 12), 1)

    layout = prefixed(etree.parse(str(code_layout)).getroot(), "gh_")
    for key in ("axes_9", "axes_10"):
        fresh = layout.find(f".//{SVG}g[@id='gh_{key}']")
        fresh.set("transform", groups[key].get("transform"))
        figure.replace(groups[key], fresh)
        groups[key] = fresh
    for defs in layout.findall(f"{SVG}defs"):
        root.append(defs)

    for axes, (title, x) in zip(("axes_9", "axes_10"), TITLES):
        found = [g for g in groups[axes].findall(f"{SVG}g") if text_of(g) == title]
        assert len(found) == 1, (axes, title)
        groups[axes].remove(found[0])
        label = etree.SubElement(figure, f"{SVG}text", x=str(x), y=str(TITLE_BASELINE), style=(
            "font-weight:700;font-size:8.4px;font-family:Arial,Helvetica,'Liberation Sans',sans-serif;"
            "text-anchor:middle;fill:#0b0b0b"))
        label.text = title
    moved = [groups["axes_9"], groups["axes_10"]]
    moved += [g for k, g in groups.items() if k.startswith("text_") and text_of(g) in ("G", "H")]
    assert len(moved) == 4
    for g in moved:
        g.set("transform", f"translate(0,{G_H_SHIFT}) {g.get('transform', '')}".strip())

    e_key = groups["axes_7"].find(f"{SVG}g[@id='legend_5']")
    e_key.set("transform", f"translate({E_KEY_SHIFT},0) {e_key.get('transform')}")

    new = prefixed(etree.parse(str(BARS)).getroot(), "ij_")
    block = etree.SubElement(figure, f"{SVG}g", id="panels_I_J_bars", transform=f"translate({X_PT},{Y_PT})")
    for child in new:
        if child.tag in (f"{SVG}defs", f"{SVG}g"):
            block.append(child)
    tree.write(str(out), xml_declaration=True, encoding="utf-8")
    subprocess.run(["inkscape", str(out), "--export-type=png", f"--export-filename={out.with_suffix('.png')}",
                    "--export-width=2600", "--export-height=3100", "--export-background=#ffffff"],
                   check=True, capture_output=True)
    print("wrote", out, "and its PNG")


if __name__ == "__main__":
    main()
