"""Build the interactive pole map and summary table for the book's main page.

Reads ``data/Laurentia_poles.csv`` (the compiled Laurentia paleomagnetic poles)
and writes two artifacts:

1. ``pole_map.ipynb`` -- a notebook holding the interactive Folium map. One
   marker per pole at the present-day sampling locality (SLAT/SLONG), colored by
   nominal age. Hovering shows the unit, age, and grade; the click popup adds the
   site and pole positions, A95, the paleolatitude of Duluth implied by the pole,
   the reference, and a link to that pole's assessment notebook. The map is
   committed as a code-cell output (input hidden) so it renders with notebook
   execution disabled. This page is the landing page's prominent "interactive
   pole map" link.
2. The Markdown pole-compilation table written into the ``index.md`` landing
   page between the ``<!-- POLE_TABLE_START -->`` / ``<!-- POLE_TABLE_END -->``
   sentinels (prose around it is hand-editable). Each unit links to its notebook
   by MyST source path, which resolves to the correct slug + deploy base.

A self-contained ``_static/Laurentia_pole_map.html`` is also written for quick
standalone preview of the map.

MyST does not serve a static ``_static/*.html`` via the ``{iframe}`` directive;
interactive Folium therefore has to come through a notebook output, which MyST
extracts to a served file and embeds as an ``<iframe src=...>``. Because that
served file's path depth is not fixed, the map's notebook links are
root-absolute and carry the deploy base path (``BASE_URL``).

Notebook links
--------------
Every pole links to a notebook in ``pole_notebooks/`` following the project's
``<age>_<Name>`` naming convention. Notebooks that already exist (possibly with
an author-chosen age prefix that differs from the compilation's nominal age, or
that bundle several polarity zones into one notebook) are mapped explicitly in
``EXISTING_NOTEBOOKS``; every other pole gets a convention-derived stem so the
link resolves once that notebook is built with the matching filename.

The map links use the MyST page *slug* (raw HTML cannot use MyST's internal
link resolver), reproduced by ``myst_slug`` from the observed MyST behavior:
filenames are lower-cased, non-alphanumeric runs become ``-``, and a leading
ordering number of three digits or fewer is stripped (e.g. ``780_Gunbarrel`` ->
``gunbarrel``) while a four-digit age is kept (``1045_Upper_Freda`` ->
``1045-upper-freda``). Links are relative to the map's own URL with
``target="_top"`` so they work both locally and under the GitHub Pages base
path.

Run with the project environment (carries pandas + folium), e.g.::

    python scripts/build_pole_map.py
"""

import json
import os
import re

import folium
import numpy as np
import pandas as pd

# --- paths -----------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CSV = os.path.join(ROOT, "data", "Laurentia_poles.csv")
MAP_HTML = os.path.join(ROOT, "_static", "Laurentia_pole_map.html")
INDEX_MD = os.path.join(ROOT, "index.md")
POLE_MAP_IPYNB = os.path.join(ROOT, "pole_map.ipynb")
PROVINCE_GEOJSON = os.path.join(ROOT, "data", "geologic_provinces",
                                "Whitmeyer2007_provinces.geojson")

TABLE_START = "<!-- POLE_TABLE_START -->"
TABLE_END = "<!-- POLE_TABLE_END -->"

# Deployed site base path. The map is raw HTML (it cannot use MyST's link
# resolver), so notebook links inside it are root-absolute and must carry the
# deploy base. The GitHub Pages workflow builds with BASE_URL=/<repo-name>;
# default to that and allow an env override (e.g. BASE_URL='' for a root-served
# local preview).
BASE_URL = os.environ.get("BASE_URL", "/2026_Laurentia_Precambrian_poles")

# Reference locality for the reported paleolatitude (Duluth, Minnesota).
DULUTH_LAT = 46.7867
DULUTH_LON = 267.8995  # 0-360 degrees east (-92.1005 deg)

# --- pole -> existing notebook stem ----------------------------------------
# Poles whose assessment notebook already exists in pole_notebooks/. The stem
# is the filename without extension. Several Osler / Mamainse polarity zones
# share one notebook.
EXISTING_NOTEBOOKS = {
    "Sept-Iles Layered Intrusion": "565_Sept_Iles",
    "Catoctin Basalts": "572_Catoctin_Basalts",
    "Callander Alkaline Complex": "575_Callander",
    "Baie des Moutons complex A": "583_Baie_des_Moutons",
    "Baie des Moutons complex B": "583_Baie_des_Moutons",
    "Franklin event grand mean": "719_Franklin_LIP",
    "Gunbarrel LIP": "780_Gunbarrel",
    "Adirondack metamorphic anorthosite": "887_Adirondack",
    "Jacobsville Formation": "990_Jacobsville",
    "Freda Sandstone": "1045_Upper_Freda",
    "Nonesuch Shale": "1078_Nonesuch",
    "Cardenas Basalts and Intrusions": "1082_Cardenas",
    "Michipicoten Island Formation": "1084_Michipicoten_Island_Formation",
    "Lake Shore Traps": "1086_Lake_Shore_Traps",
    "Central Arizona diabases -N": "1088_Central_Arizona_intrusions",
    "Schroeder Lutsen Basalts": "1090_Schroeder_Lutsen_Basalts",
    "Portage Lake Volcanics": "1092_Portage_Lake_Volcanics",
    "North Shore lavas -N": "1095_North_Shore_Volcanic_Group",
    "Upper Osler volcanics -R": "1108_Osler_Volcanic_Group",
    "Middle Osler volcanics -R": "1108_Osler_Volcanic_Group",
    "Lower Osler volcanics -R": "1108_Osler_Volcanic_Group",
    "Uppermost Mamainse Point volcanics -N": "1109_Mamainse_Point_Volcanics",
    "Mamainse Point volcanics -C (lower N, upper R)": "1109_Mamainse_Point_Volcanics",
    "Lower Mamainse Point volcanics -R2": "1109_Mamainse_Point_Volcanics",
    "Lowermost Mamainse Point volcanics -R1": "1109_Mamainse_Point_Volcanics",
    "Cleaver Dykes": "1740_Cleaver_Dikes",
    "NE trending ECMB Diabase Dykes": "1779_East_Central_Minnesota_Batholith",
    # Notebooks added with short custom stems (do not match the <age>_<Name>
    # convention), mapped explicitly so the table/map links resolve.
    "Kwagunt Formation": "757_Chuar_Group",  # rebuilt as the combined Chuar Group pole
    "Uinta Mountain Group": "754_Uinta",
    "Abitibi Dykes": "1141_Abitibi",
    "NW Ontario Lamprophyre Dykes and Abitibi Dykes": "1144_Lamprophyre",
    "Sudbury Dykes Combined": "1237_Sudbury",
    "Mackenzie dykes grand mean": "1267_Mackenzie",
    "Nain Anorthosite": "1305_Nain",
    "Midsommersoe Dolerite": "1382_Midsommersoe",
    "Victoria Fjord dolerite dykes": "1382_Victoria",
    "Zig-Zag Dal Basalts": "1382_Zigzag",
    "MEAN Rocky Mountain intrusions": "1430_Rocky",
    "Michikamau Intrusion Combined": "1460_Michikamau",
    "St.Francois Mountains Acidic Rocks": "1476_Francois",
    "Melville Bugt diabase dykes": "1633_Melville",
}


def convention_stem(rockname, age):
    """Convention-derived notebook stem ``<age>_<SanitizedName>`` for a pole.

    Used for poles that do not yet have a notebook; the eventual notebook
    should use this filename so its link (and map slug) resolves.

    Args:
        rockname (str): The compilation ``ROCKNAME``.
        age (int): The nominal age in Ma.

    Returns:
        str: The notebook stem (filename without extension).
    """
    name = re.sub(r"[^0-9A-Za-z]+", "_", str(rockname)).strip("_")
    return f"{age}_{name}"


def _name_tokens(text):
    """Lowercase alphanumeric tokens of a string."""
    return {t for t in re.split(r"[^a-z0-9]+", str(text).lower()) if t}


def build_stem_map(df):
    """Map each pole's ROCKNAME to its notebook stem.

    Notebook filenames use short, age-prefixed abbreviations of the unit name
    (e.g. ``2218_Senneterre``). Resolution order per pole:

    1. an explicit entry in :data:`EXISTING_NOTEBOOKS` (multi-pole notebooks and
       names whose tokens don't match, e.g. ``Zig-Zag`` -> ``Zigzag``);
    2. a notebook file on disk whose abbreviated name tokens are all contained
       in the ROCKNAME, choosing the closest age when several poles qualify;
    3. the ``<age>_<Name>`` convention stem (for poles whose notebook is not yet
       built).

    Args:
        df (pd.DataFrame): The loaded pole table (needs ROCKNAME, nominal age).

    Returns:
        dict: ROCKNAME -> notebook stem.
    """
    nb_dir = os.path.join(ROOT, "pole_notebooks")
    parsed = []
    for stem in sorted(f[:-6] for f in os.listdir(nb_dir)
                       if f.endswith(".ipynb")):
        mo = re.match(r"^(\d+)[_-](.*)$", stem)
        age = int(mo.group(1)) if mo else None
        parsed.append((stem, age, _name_tokens(mo.group(2) if mo else stem)))

    poles = list(zip(df["ROCKNAME"], df["nominal age"].astype(int)))
    pole_tokens = {rn: _name_tokens(rn) for rn, _ in poles}

    mapping = {rn: EXISTING_NOTEBOOKS[rn] for rn, _ in poles
               if rn in EXISTING_NOTEBOOKS}
    claimed = set(mapping)
    used = set(EXISTING_NOTEBOOKS.values())

    # Match more-specific (more-token) notebooks first to resolve ambiguity.
    for stem, age, ntoks in sorted(parsed, key=lambda x: -len(x[2])):
        if stem in used or not ntoks:
            continue
        cands = [(rn, page) for rn, page in poles
                 if rn not in claimed and ntoks <= pole_tokens[rn]]
        if not cands:
            continue
        rn = min(cands,
                 key=lambda c: abs(c[1] - (age if age is not None else c[1])))[0]
        mapping[rn] = stem
        claimed.add(rn)
        used.add(stem)

    for rn, age in poles:
        mapping.setdefault(rn, convention_stem(rn, age))
    return mapping


def myst_slug(stem):
    """Reproduce the MyST page slug for a ``pole_notebooks/<stem>.ipynb`` file.

    MyST lower-cases the filename, turns non-alphanumeric runs into ``-``, and
    strips a leading ordering number of three digits or fewer (a four-digit
    age is retained). Verified against the built ``myst.xref.json``.

    Args:
        stem (str): Notebook filename without extension.

    Returns:
        str: The page slug used in the deployed URL.
    """
    s = stem
    mo = re.match(r"^(\d{1,3})[_\-\s]+(.*)$", s)
    if mo:
        s = mo.group(2)
    s = re.sub(r"[^a-z0-9]+", "-", s.lower())
    return s.strip("-")


def paleolatitude(site_lat, site_lon, pole_lat, pole_lon):
    """Paleolatitude of a locality implied by a paleomagnetic pole (GAD).

    The angular distance between the locality and the pole is the paleo
    co-latitude, so the paleolatitude is ``90 - distance``.

    Args:
        site_lat, site_lon (float): Locality latitude and longitude (deg,
            longitude in 0-360 deg east).
        pole_lat, pole_lon (float): Pole latitude and longitude (deg).

    Returns:
        float: Paleolatitude of the locality in degrees.
    """
    slat, slon = np.radians(site_lat), np.radians(site_lon)
    plat, plon = np.radians(pole_lat), np.radians(pole_lon)
    cos_p = (np.sin(slat) * np.sin(plat)
             + np.cos(slat) * np.cos(plat) * np.cos(plon - slon))
    return 90.0 - np.degrees(np.arccos(np.clip(cos_p, -1.0, 1.0)))


def short_reference(authors, year):
    """Compact 'first-author et al. (year)' reference string.

    The source CSV stores a few author fields with a corrupted character
    (U+FFFD) where 'Luleå' was lost; this is repaired to the known Nordic
    "Luleå Working Group Mean" attribution.
    """
    authors = "" if pd.isna(authors) else str(authors).strip()
    if "LULE�" in authors.upper() or "� WORKING GROUP" in authors:
        return f"Luleå Working Group Mean ({fmt_year(year)})"
    surname = authors.split(",")[0].strip() if authors else "—"
    etal = " et al." if authors.count(",") > 1 else ""
    return f"{surname}{etal} ({fmt_year(year)})"


def fmt_year(year):
    """Format a publication year, or '—' when missing."""
    try:
        return f"{int(year)}"
    except (TypeError, ValueError):
        return "—"


def fmt(value, decimals=1):
    """Format a number to fixed decimals, or '–' when missing."""
    try:
        x = float(value)
    except (TypeError, ValueError):
        return "–"
    return "–" if np.isnan(x) else f"{x:.{decimals}f}"


def grade_label(grade):
    """Reliability grade (A/B) as a string, or '–' when missing."""
    return "–" if pd.isna(grade) else str(grade).strip()


def to_pm180(lon):
    """Convert a 0-360 deg longitude to the -180/180 convention for Folium."""
    return ((float(lon) + 180.0) % 360.0) - 180.0


def age_colormap(df):
    """Continuous viridis color scale over the compilation's age range."""
    import branca.colormap as bcm
    import matplotlib.cm as mcm
    import matplotlib.colors as mcolors

    colors = [mcolors.to_hex(mcm.viridis(i / 255)) for i in range(0, 256, 16)]
    cmap = bcm.LinearColormap(
        colors,
        vmin=float(df["nominal age"].min()),
        vmax=float(df["nominal age"].max()),
        caption="Nominal age (Ma)",
    )
    return cmap


def load_provinces():
    """Load the basement-province GeoJSON, or None if it has not been built."""
    if not os.path.exists(PROVINCE_GEOJSON):
        return None
    with open(PROVINCE_GEOJSON) as fh:
        return json.load(fh)


def _province_legend_html(provinces):
    """Fixed-position HTML legend for province categories and pole markers."""
    rows = []
    for feat in provinces["features"]:
        p = feat["properties"]
        rows.append(
            f'<div><span style="display:inline-block;width:13px;height:13px;'
            f'background:{p["color"]};border:1px solid #777;margin-right:6px;'
            f'vertical-align:middle;"></span>{p["category"]}</div>')
    poles = (
        '<div style="margin-top:5px;border-top:1px solid #ccc;padding-top:4px;">'
        '<span style="display:inline-block;width:12px;height:12px;'
        'border-radius:50%;background:#555;border:1px solid #000;'
        'margin-right:6px;vertical-align:middle;"></span>Grade A pole<br>'
        '<span style="display:inline-block;width:11px;height:11px;'
        'background:#555;border:1px solid #000;margin-right:6px;'
        'vertical-align:middle;"></span>Grade B pole</div>')
    return (
        '<div style="position: fixed; bottom: 22px; left: 12px; z-index: 1000;'
        ' background: rgba(255,255,255,0.9); padding: 8px 11px;'
        ' border: 1px solid #999; border-radius: 4px; font-size: 12px;'
        ' line-height: 1.5; max-width: 230px;">'
        '<div style="font-weight:600;margin-bottom:3px;">Basement provinces'
        '</div>' + "".join(rows) + poles +
        '<div style="margin-top:5px;font-size:10px;color:#666;">Provinces after'
        ' Whitmeyer &amp; Karlstrom (2007)</div></div>')


def _pole_marker(row, color, slug):
    """Build a clickable pole marker: circle for Grade A, square for Grade B."""
    slat, slon = float(row["SLAT"]), to_pm180(row["SLONG"])
    grade = grade_label(row["Grade"])
    href = f"{BASE_URL}/{slug}"
    tooltip = (f'<div style="font-size:17px;line-height:1.4;">'
               f"<b>{row['ROCKNAME']}</b><br>"
               f"~{int(row['nominal age'])} Ma · Grade {grade}</div>")
    popup_html = (
        f'<div style="font-size:16px;line-height:1.5;">'
        f"<b>{row['ROCKNAME']}</b><br>"
        f"<i>{row['Terrane']}</i><br>"
        f"~{int(row['nominal age'])} Ma &middot; Grade {grade}<br>"
        f"Site: {fmt(row['SLAT'])}&deg;N, {fmt(slon)}&deg;E<br>"
        f"Pole: {fmt(row['PLAT'])}&deg;N, {fmt(row['PLONG'])}&deg;E "
        f"(A95 {fmt(row['A95'])}&deg;)<br>"
        f"Duluth paleolat: {fmt(row['_duluth_paleolat'])}&deg;<br>"
        f"{short_reference(row['POLE AUTHORS'], row['YEAR'])}<br>"
        f"<a href='{href}' target='_blank' rel='noopener'>"
        f"Open notebook &rarr;</a></div>")
    popup = folium.Popup(popup_html, max_width=320)

    if grade == "A":
        return folium.CircleMarker(
            location=[slat, slon], radius=6, color="black", weight=1,
            fill=True, fill_color=color, fill_opacity=0.9,
            tooltip=tooltip, popup=popup)
    square = (f'<div style="width:12px;height:12px;background:{color};'
              f'border:1.5px solid black;transform:translate(-50%,-50%);">'
              f'</div>')
    return folium.Marker(
        [slat, slon],
        icon=folium.DivIcon(icon_size=(0, 0), icon_anchor=(0, 0), html=square),
        tooltip=tooltip, popup=popup)


def build_map(df):
    """Build the enriched interactive Folium map.

    Layers (toggleable via the layer control): basement geologic provinces
    (Whitmeyer & Karlstrom, 2007) colored by age, their province/orogen name
    labels, and the paleomagnetic-pole sampling localities. Each pole marker is
    placed at SLAT/SLONG, colored by nominal age, drawn as a circle (Grade A) or
    square (Grade B), and is clickable for a popup with the pole position, A95,
    Duluth paleolatitude, reference, and a link to its assessment notebook.
    """
    cmap = age_colormap(df)
    # prefer_canvas renders the vector layers (provinces, circle markers) on a
    # canvas, which paints reliably inside the notebook output iframe (the SVG
    # renderer can fail to paint there even though the data is present).
    m = folium.Map(location=[52, -90], zoom_start=4, tiles=None,
                   world_copy_jump=True, prefer_canvas=True)
    folium.TileLayer("CartoDB positron", name="Light basemap",
                     control=True).add_to(m)
    folium.TileLayer("OpenStreetMap", name="OpenStreetMap",
                     show=False).add_to(m)

    # Enlarge the layer-control ("map selection") and colorbar text. The layer
    # control class is created by Leaflet at runtime; the colorbar is a branca
    # SVG whose id begins with "color_map_".
    m.get_root().header.add_child(folium.Element(
        "<style>"
        ".leaflet-control-layers, .leaflet-control-layers label,"
        " .leaflet-control-layers span { font-size: 15px; line-height: 1.6; }"
        ' svg[id^="color_map_"] text { font-size: 15px; }'
        "</style>"))

    provinces = load_provinces()
    if provinces is not None:
        prov_group = folium.FeatureGroup(name="Geologic provinces", show=True)
        folium.GeoJson(
            provinces,
            style_function=lambda feat: {
                "fillColor": feat["properties"]["color"],
                "color": "#777", "weight": 0.4, "fillOpacity": 0.55},
        ).add_to(prov_group)
        prov_group.add_to(m)

    pole_group = folium.FeatureGroup(name="Paleomagnetic poles", show=True)
    for _, row in df.iterrows():
        _pole_marker(row, cmap(float(row["nominal age"])),
                     myst_slug(row["_stem"])).add_to(pole_group)
    pole_group.add_to(m)

    cmap.add_to(m)
    if provinces is not None:
        m.get_root().html.add_child(
            folium.Element(_province_legend_html(provinces)))
    folium.LayerControl(collapsed=False).add_to(m)
    # Frame Laurentia (Duluth-centered extent of the static compilation map).
    m.fit_bounds([[27, -119], [78, -57]])
    return m


def build_table(df):
    """Build the Markdown pole-compilation table for the landing page.

    Each unit links to its notebook by MyST source path
    (``pole_notebooks/<stem>.ipynb``), which MyST resolves to the correct
    slug + deploy base for any notebook in the toc. "Rating" is the Nordic
    grade (A/B).

    Args:
        df (pd.DataFrame): Output of :func:`load_poles`.

    Returns:
        str: A GitHub-flavored Markdown table.
    """
    header = ("| Terrane | Unit | Age (Ma) | Rating | Site lon | Site lat "
              "| Plon | Plat | A95 | Duluth paleolat | Pole reference |")
    sep = "|" + "|".join(["---"] * 11) + "|"
    lines = [header, sep]
    for _, row in df.iterrows():
        link = f"[{row['ROCKNAME']}](pole_notebooks/{row['_stem']}.ipynb)"
        lines.append(
            f"| {row['Terrane']} | {link} | {int(row['nominal age'])} "
            f"| {grade_label(row['Grade'])} "
            f"| {fmt(row['SLONG'])} | {fmt(row['SLAT'])} "
            f"| {fmt(row['PLONG'])} | {fmt(row['PLAT'])} | {fmt(row['A95'])} "
            f"| {fmt(row['_duluth_paleolat'])} "
            f"| {short_reference(row['POLE AUTHORS'], row['YEAR'])} |"
        )
    return "\n".join(lines)


def load_poles():
    """Load and prepare the compilation table for the map and table.

    Returns:
        pd.DataFrame: ``data/Laurentia_poles.csv`` with usable rows only and the
        added helper columns ``_stem`` (notebook filename stem) and
        ``_duluth_paleolat`` (paleolatitude of Duluth implied by the pole).
    """
    df = pd.read_csv(CSV)
    df["nominal age"] = pd.to_numeric(df["nominal age"], errors="coerce")
    df = df.dropna(subset=["SLAT", "SLONG", "PLAT", "PLONG",
                           "nominal age"]).copy()
    df["nominal age"] = df["nominal age"].astype(int)
    df["_stem"] = df["ROCKNAME"].map(build_stem_map(df))
    df["_duluth_paleolat"] = paleolatitude(
        DULUTH_LAT, DULUTH_LON, df["PLAT"].astype(float),
        df["PLONG"].astype(float))
    return df


# Markdown for the interactive-map notebook's lead cell.
MAP_NOTEBOOK_MD = """\
# Interactive pole map

An interactive map of the compilation over the basement geologic provinces of \
Laurentia (after Whitmeyer & Karlstrom, 2007). Each marker is the present-day \
sampling locality of a pole, colored by nominal age and drawn as a circle \
(Grade A) or square (Grade B). Hover for the unit, age, and grade; click a \
marker for the site and pole positions, A95, the paleolatitude it implies for \
Duluth, Minnesota, the reference, and a link to that pole's assessment \
notebook. Use the layer control (top right) to toggle the provinces and poles, \
and zoom/pan freely. Built from `data/Laurentia_poles.csv`."""

# Source for the map code cell. Input is removed in the rendered book; this is
# what re-executing the cell against the current CSV renders inline.
MAP_CELL_SOURCE = """\
import sys
sys.path.insert(0, "scripts")
from build_pole_map import load_poles, build_map

# Display the Folium map object directly so it renders live in the notebook;
# its root-absolute notebook links still resolve once deployed.
build_map(load_poles())"""


def write_pole_map_notebook(df):
    """Write ``pole_map.ipynb`` with the interactive map as a baked code cell.

    The map code cell carries a committed HTML output (the Folium map's
    ``_repr_html_``) so the page renders the interactive map with notebook
    execution disabled; re-running the cell (or this script) rebuilds it from
    the current ``data/Laurentia_poles.csv``.
    """
    import nbformat as nbf

    map_cell = nbf.v4.new_code_cell(source=MAP_CELL_SOURCE)
    map_cell.metadata = {"tags": ["remove-input"]}
    map_cell.execution_count = 1
    map_cell.outputs = [nbf.v4.new_output(
        "execute_result", data={"text/html": build_map(df)._repr_html_()},
        metadata={}, execution_count=1)]

    nb = nbf.v4.new_notebook()
    nb.cells = [nbf.v4.new_markdown_cell(MAP_NOTEBOOK_MD), map_cell]
    nb.metadata["kernelspec"] = {
        "display_name": "Python 3", "language": "python", "name": "python3"}
    with open(POLE_MAP_IPYNB, "w", encoding="utf-8") as fh:
        nbf.write(nb, fh)


def write_table_into_index(table_md):
    """Replace the table block in index.md between the sentinel comments."""
    with open(INDEX_MD, encoding="utf-8") as fh:
        text = fh.read()
    block = f"{TABLE_START}\n{table_md}\n{TABLE_END}"
    if TABLE_START not in text or TABLE_END not in text:
        raise SystemExit(
            f"Sentinels {TABLE_START}/{TABLE_END} not found in index.md; "
            "add them where the table should go.")
    text = re.sub(re.escape(TABLE_START) + r".*?" + re.escape(TABLE_END),
                  block, text, flags=re.DOTALL)
    with open(INDEX_MD, "w", encoding="utf-8") as fh:
        fh.write(text)


def report_link_coverage(df):
    """Report how pole links map onto the notebooks that actually exist.

    Flags two kinds of drift so naming mismatches surface at build time:
    poles whose linked stem has no file on disk, and notebooks present in
    ``pole_notebooks/`` that no pole row links to.
    """
    nb_dir = os.path.join(ROOT, "pole_notebooks")
    on_disk = {f[:-6] for f in os.listdir(nb_dir) if f.endswith(".ipynb")}
    linked = {stem: rn for stem, rn in zip(df["_stem"], df["ROCKNAME"])}

    resolved = sum(stem in on_disk for stem in linked)
    print(f"{resolved}/{len(linked)} pole links resolve to an existing notebook.")

    missing = [(rn, stem) for stem, rn in linked.items() if stem not in on_disk]
    explicit_missing = [(rn, st) for rn, st in missing
                        if rn in EXISTING_NOTEBOOKS]
    if explicit_missing:
        print("  WARNING: explicitly-mapped notebooks not found on disk:")
        for rn, st in explicit_missing:
            print(f"    {rn!r} -> pole_notebooks/{st}.ipynb")

    orphans = sorted(on_disk - set(linked))
    if orphans:
        print("  Notebooks on disk that no pole row links to (check naming):")
        for st in orphans:
            print(f"    pole_notebooks/{st}.ipynb")


def main():
    df = load_poles()

    os.makedirs(os.path.dirname(MAP_HTML), exist_ok=True)
    build_map(df).save(MAP_HTML)          # standalone copy for direct preview
    write_pole_map_notebook(df)
    write_table_into_index(build_table(df))

    print(f"Wrote {os.path.relpath(MAP_HTML, ROOT)} (preview, {len(df)} poles)")
    print(f"Wrote {os.path.relpath(POLE_MAP_IPYNB, ROOT)} (interactive map)")
    print(f"Updated table block in {os.path.relpath(INDEX_MD, ROOT)} "
          f"({len(df)} rows)")
    print(f"Map notebook links use BASE_URL={BASE_URL!r}")
    report_link_coverage(df)


if __name__ == "__main__":
    main()
