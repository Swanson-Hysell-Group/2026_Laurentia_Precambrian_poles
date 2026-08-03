"""Build the interactive pole map and summary table for the book's main page.

Poles come from two sources, disjoint in age so their union introduces no
duplicates. Poles at or below ``PUBLISH_AGE_MAX`` (1779 Ma) are read from
``data/nordic_summaries/nordic_summaries_combined.csv``, the per-notebook
summaries that ``data/nordic_summaries/combine_nordic_summaries.py`` assembles
from the assessment notebooks in ``pole_notebooks/``; these are the poles
recreated at the site level as part of this project, so the map and table track
the current notebook values. Older poles are carried over from
``data/Laurentia_poles.csv``, the compilation assembled at earlier Nordic
Paleomagnetic Workshops, and are listed without a notebook link. From that
combined set the script writes two artifacts:

1. ``pole_map.ipynb`` -- a notebook holding the interactive Folium map. One
   marker per pole at the present-day sampling locality (SLAT/SLONG), colored by
   nominal age. Hovering shows the unit, age, and grade; the click popup adds the
   site and pole positions, A95, the paleolatitude of Duluth implied by the pole
   (where the pole constrains it -- see :data:`DULUTH_PALEOLAT_TERRANES`), the
   reference, and a link to that pole's assessment notebook. The map is
   committed as a code-cell output (input hidden) so it renders with notebook
   execution disabled. This page is the landing page's prominent "interactive
   pole map" link.
2. The Markdown pole-compilation table written into the ``compilation.md`` page
   between the ``<!-- POLE_TABLE_START -->`` / ``<!-- POLE_TABLE_END -->``
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

import glob
import json
import os
import re

import folium
import numpy as np
import pandas as pd
import pmagpy.pmag as pmag

# --- paths -----------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CSV = os.path.join(ROOT, "data", "Laurentia_poles.csv")
# Per-notebook Nordic summaries (the latest recreated poles) and their combined
# table. The map and compilation table are built from these for the published
# (<= PUBLISH_AGE_MAX) poles, so they reflect the current notebook values; older
# poles are carried from the legacy compilation CSV above.
SUMMARY_DIR = os.path.join(ROOT, "data", "nordic_summaries")
COMBINED_CSV = os.path.join(SUMMARY_DIR, "nordic_summaries_combined.csv")
# Inclination-shallowing-corrected (Kent) mean poles for the sedimentary units,
# assembled by data/nordic_summaries/build_kent_poles.py. These positions
# replace the uncorrected ones for those units wherever a pole is plotted.
KENT_CSV = os.path.join(SUMMARY_DIR, "kent_poles_combined.csv")
MAP_HTML = os.path.join(ROOT, "_static", "Laurentia_pole_map.html")
COMPILATION_MD = os.path.join(ROOT, "compilation.md")
POLE_MAP_IPYNB = os.path.join(ROOT, "pole_map.ipynb")
PROVINCE_GEOJSON = os.path.join(ROOT, "data", "geologic_provinces",
                                "Whitmeyer2007_provinces.geojson")

TABLE_START = "<!-- POLE_TABLE_START -->"
TABLE_END = "<!-- POLE_TABLE_END -->"
KENT_TABLE_START = "<!-- KENT_TABLE_START -->"
KENT_TABLE_END = "<!-- KENT_TABLE_END -->"

# Deployed site base path. The map is raw HTML (it cannot use MyST's link
# resolver), so notebook links inside it are root-absolute and must carry the
# deploy base. The GitHub Pages workflow builds with BASE_URL=/<repo-name>;
# default to that and allow an env override (e.g. BASE_URL='' for a root-served
# local preview).
BASE_URL = os.environ.get("BASE_URL", "/2026_Laurentia_Precambrian_poles")

# Only assessment notebooks with a stem prefix at or below this age (Ma) are
# published in the book. Poles whose notebook stem prefix is larger are still
# shown on the map and listed in the table, but without a notebook link.
PUBLISH_AGE_MAX = 1779

# Notebooks that exist but are deliberately not published yet, so the
# compilation table and map list the pole without linking to a page that does
# not exist in the book. Currently the Wharton Group notebook, which
# transcribes site-level data from the unpublished Raub et al. (2026)
# manuscript; the pole summary itself is public, the site table is not. Remove
# the entry and add the notebook to myst.yml once that paper is out.
WITHHELD_NOTEBOOKS = frozenset({"1756_Wharton_Group"})

# Reference locality for the reported paleolatitude (Duluth, Minnesota).
DULUTH_LAT = 46.7867
DULUTH_LON = 267.8995  # 0-360 degrees east (-92.1005 deg)

# Terranes that have rifted from Laurentia, rotated back into the Laurentia
# reference frame before the Duluth paleolatitude is computed so that column is
# comparable across the whole table. Mirrors pole_tools.TERRANE_EULER_POLES
# (Greenland: Roest & Srivastava, 1989; Scotland: Torsvik & Cocks, 2017). The
# tabulated Plat/Plon stay present-day; only the paleolatitude is rotated.
TERRANE_EULER_POLES = {
    "Laurentia-Greenland": [67.5, -118.5, -13.8],
    "Laurentia-Greenland-Nain": [67.5, -118.5, -13.8],
    "Laurentia-Scotland": [78.6, 161.9, -31.0],
    "Laurentia-Svalbard": [-81.0, 125.0, 68.0],
}

# Duluth sits on the western Superior craton. Prior to the assembly of
# Laurentia, a pole from a block that had not yet joined places that block, not
# Duluth, so the Duluth paleolatitude is only reported for poles older than
# PUBLISH_AGE_MAX whose terrane is the Superior craton itself or the
# Trans-Hudson orogen, which welded to it. Poles from the Slave, Rae, Wyoming,
# Nain, and eastern Superior blocks over that interval are left blank rather
# than implying a Duluth position that the pole does not constrain. Mirrors the
# terranes that pole_tools.get_Laurentia_poles treats as being in the Laurentia
# reference frame.
DULUTH_PALEOLAT_TERRANES = frozenset({
    "Laurentia-Superior",
    "Laurentia-Superior(West)",
    "Laurentia-Trans-Hudson orogen",
})

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
    "Freda Sandstone": "1075_Lower_Freda",  # legacy compilation pole = Henry et al. (1977) lower Freda
    "Nonesuch Formation": "1078_Nonesuch",
    "Cardenas Basalts and Intrusions": "1082_Cardenas",
    "Michipicoten Island Formation": "1084_Michipicoten_Island_Formation",
    "Lake Shore Traps": "1086_Lake_Shore_Traps",
    "Central Arizona diabases -N": "1098_Central_Arizona_intrusions",
    "Schroeder Lutsen Basalts": "1090_Schroeder_Lutsen_Basalts",
    "Portage Lake Volcanics": "1092_Portage_Lake_Volcanics",
    "North Shore lavas -N": "1095_North_Shore_Volcanic_Group",
    "Chengwatana Volcanics": "1096_Chengwatana",
    "MEAN Nipigon sills and lavas": "1108_Nipigon",
    "North Qoroq Intr.": "1275_North_Qoroq",
    "Kungnat Ring Dyke": "1275_Kungnat",
    "South Qoroq Intr.": "1163_South_Qoroq",
    "NE-SW Trending Dyke Swarm": "1160_NE-SW",
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
    "Kwagunt Formation": "755_Chuar_Group",  # rebuilt as the combined Chuar Group pole
    "Uinta Mountain Group": "759_Uinta_Mountain_Group",
    "Abitibi Dykes": "1141_Abitibi",
    "NW Ontario Lamprophyre Dykes and Abitibi Dykes": "1144_Lamprophyre",
    "Sudbury Dykes Combined": "1235_Sudbury",
    "Mackenzie dykes grand mean": "1267_Mackenzie",
    "Nain Anorthosite": "1305_Nain",
    "Midsommersoe Dolerite": "1382_Midsommersoe",
    "Victoria Fjord dolerite dykes": "1382_Victoria",
    "Zig-Zag Dal Basalts": "1382_Zigzag",
    "MEAN Rocky Mountain intrusions": "1430_Rocky",
    "Michikamau Intrusion Combined": "1469_Michikamau",
    "St.Francois Mountains Acidic Rocks": "1466_Francois",
    "Melville Bugt diabase dykes": "1630_Melville",
    "Elbow Creek dikes": "2480_Elbow_creek",  # not published (>1779), mapped so the stem resolves
}

# Per-pole summary CSVs whose filename stem is not itself a notebook filename.
# Several notebooks emit more than one summary (the Mamainse and Osler polarity
# zones, the two Baie des Moutons complexes), and a few summaries are named
# slightly differently from their notebook; map each summary stem to the
# notebook that should be linked.
SUMMARY_STEM_TO_NOTEBOOK = {
    "583_Baie_des_Moutons_A": "583_Baie_des_Moutons",
    "583_Baie_des_Moutons_B": "583_Baie_des_Moutons",
    "1094_Mamainse_upper_N": "1109_Mamainse_Point_Volcanics",
    "1100_Mamainse_Flour_Bay": "1109_Mamainse_Point_Volcanics",
    "1105_Mamainse_lower_R2": "1109_Mamainse_Point_Volcanics",
    "1109_Mamainse_lower_R1": "1109_Mamainse_Point_Volcanics",
    "1105_Osler_reverse_upper": "1108_Osler_Volcanic_Group",
    "1107_Osler_reverse_middle": "1108_Osler_Volcanic_Group",
    "1108_Osler_reverse_lower": "1108_Osler_Volcanic_Group",
    "1382_Midsommersoe_Dolerites": "1382_Midsommersoe",
    "1382_Victoria_Fjord": "1382_Victoria",
    "1382_ZigZag_Dal_Basalt": "1382_Zigzag",
    "1430_Rocky_Mountain_intrusions": "1430_Rocky",
    "1466_St_Francois_Mountains": "1466_Francois",
    "1592_WesternChannelDiabase": "1592_Western_Channel",
    "1630_Melville_Bugt": "1630_Melville",
}


def summary_stem_by_rockname():
    """Map each summary ROCKNAME to the notebook stem that should be linked.

    Reads the per-pole Nordic summary CSVs that compose ``COMBINED_CSV`` (each
    filename is the summary stem) and resolves it to a notebook stem via
    :data:`SUMMARY_STEM_TO_NOTEBOOK`, defaulting to the summary stem itself when
    that is already a notebook filename.
    """
    out = {}
    for path in sorted(glob.glob(os.path.join(SUMMARY_DIR, "*.csv"))):
        stem = os.path.basename(path)[:-4]
        if not stem[:1].isdigit() or stem == "nordic_summaries_combined":
            continue
        df1 = pd.read_csv(path)
        if df1.empty or "ROCKNAME" not in df1.columns:
            continue
        rockname = str(df1.iloc[0]["ROCKNAME"]).strip()
        out[rockname] = SUMMARY_STEM_TO_NOTEBOOK.get(stem, stem)
    return out


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


def stem_age(stem):
    """Leading age prefix (Ma) of a notebook stem, or None if it has none."""
    mo = re.match(r"^(\d+)", str(stem))
    return int(mo.group(1)) if mo else None


def is_published(stem):
    """True if the notebook stem is published (prefix <= PUBLISH_AGE_MAX, not withheld)."""
    if stem in WITHHELD_NOTEBOOKS:
        return False
    age = stem_age(stem)
    return age is not None and age <= PUBLISH_AGE_MAX


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


# Kent 95% ellipse parameters carried through ``apply_kent_poles`` (prefixed
# with "_" on the DataFrame) so a map can draw the true ellipse rather than the
# equal-area circle. Zdec/Zinc is the major-axis direction and Zeta its
# semi-angle; Edec/Einc and Eta are the minor axis.
KENT_ELLIPSE_COLUMNS = ("Zdec", "Zinc", "Zeta", "Edec", "Einc", "Eta")


def kent_dict(row):
    """Repackage a row's Kent columns into an ``ipmag.plot_pole_ellipse`` dict.

    Args:
        row (pd.Series): A row carrying ``PLAT``/``PLONG`` and the ``_Zdec``…
            ``_Eta`` columns added by :func:`apply_kent_poles`.

    Returns:
        dict: ``{'dec', 'inc', 'Zdec', 'Zinc', 'Zeta', 'Edec', 'Einc', 'Eta'}``.
    """
    d = {"dec": float(row["PLONG"]), "inc": float(row["PLAT"])}
    d.update({c: float(row["_" + c]) for c in KENT_ELLIPSE_COLUMNS})
    return d


def load_kent_poles():
    """The Kent mean poles of the sedimentary units, keyed by ROCKNAME.

    Reads ``KENT_CSV``. Returns an empty mapping (with a warning) if that file
    has not been built yet, so the map and figures still render from the
    uncorrected positions alone.
    """
    if not os.path.exists(KENT_CSV):
        print(f"-W- {os.path.basename(KENT_CSV)} not found; sedimentary poles "
              "are plotted uncorrected. Run "
              "data/nordic_summaries/build_kent_poles.py")
        return {}
    kent = pd.read_csv(KENT_CSV)
    return {str(row["ROCKNAME"]).strip(): row for _, row in kent.iterrows()}


def apply_kent_poles(df, verbose=True):
    """Substitute the inclination-shallowing-corrected pole for sedimentary units.

    Detrital remanence is shallowed during compaction, so the as-measured pole
    of a sedimentary unit implies a paleolatitude that is a minimum. For every
    unit in the Kent table (see
    ``data/nordic_summaries/build_kent_poles.py``) the ``PLAT``/``PLONG``
    columns are replaced by the Kent mean pole and ``A95`` by the equal-area
    circular approximation to its 95% confidence ellipse,
    ``sqrt(zeta95 * eta95)`` — the radius consumers that can only draw a
    circular confidence should use. The uncorrected values are preserved in
    ``_plat_uncorrected`` / ``_plong_uncorrected`` / ``_a95_uncorrected`` so a
    popup or caption can still report them, and ``_kent`` flags the substituted
    rows.

    The full ellipse is carried alongside in :data:`KENT_ELLIPSE_COLUMNS`
    (``_Zdec``, ``_Zinc``, ``_Zeta``, ``_Edec``, ``_Einc``, ``_Eta``) so a map
    can draw the real Kent ellipse with ``ipmag.plot_pole_ellipse`` instead of
    the circular approximation; :func:`kent_dict` packages a row's values back
    into the dictionary that function expects.

    Args:
        df (pd.DataFrame): Poles with ``ROCKNAME``, ``PLAT``, ``PLONG``,
            ``A95``.
        verbose (bool): Print which units were substituted.

    Returns:
        pd.DataFrame: A copy with the corrected positions and the added
        ``_kent``, ``_kent_method``, ``_*_uncorrected``, and ellipse columns.
    """
    df = df.copy()
    df["_kent"] = False
    df["_kent_method"] = ""
    for col, src in (("_plat_uncorrected", "PLAT"),
                     ("_plong_uncorrected", "PLONG"),
                     ("_a95_uncorrected", "A95")):
        df[col] = pd.to_numeric(df[src], errors="coerce")
    for col in KENT_ELLIPSE_COLUMNS:
        df["_" + col] = np.nan

    kent = load_kent_poles()
    substituted = []
    for idx in df.index:
        row = kent.get(str(df.at[idx, "ROCKNAME"]).strip())
        if row is None:
            continue
        df.at[idx, "PLAT"] = float(row["PLAT"])
        df.at[idx, "PLONG"] = float(row["PLONG"])
        df.at[idx, "A95"] = float(row["A95"])
        df.at[idx, "_kent"] = True
        df.at[idx, "_kent_method"] = str(row["f method"])
        for col in KENT_ELLIPSE_COLUMNS:
            df.at[idx, "_" + col] = float(row[col])
        substituted.append(str(df.at[idx, "ROCKNAME"]))
    if verbose and substituted:
        print(f"  Kent (inclination-corrected) poles used for "
              f"{len(substituted)} sedimentary unit(s): "
              + ", ".join(sorted(substituted)))
    unused = sorted(set(kent) - set(substituted))
    if verbose and unused:
        print("  -W- Kent poles with no matching row here: "
              + ", ".join(unused))
    return df


def short_reference(authors, year):
    """Compact 'first-author et al. (year)' reference string.

    Some source author fields carry a degraded form of the Nordic "Luleå
    Working Group Mean" attribution (e.g. 'LULE? WORKING GROUP MEAN' where the
    'å' was lost to a U+FFFD or dropped); these are repaired to the proper name.
    """
    authors = "" if pd.isna(authors) else str(authors).strip()
    if "LULE" in authors.upper() and "WORKING GROUP" in authors.upper():
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
        '<br>Whitmeyer &amp; Karlstrom (2007)</div></div>')


def _pole_marker(row, color, slug):
    """Build a clickable pole marker: circle for Grade A, square for Grade B."""
    slat, slon = float(row["SLAT"]), to_pm180(row["SLONG"])
    grade = grade_label(row["Grade"])
    href = f"{BASE_URL}/{slug}"
    tooltip = (f'<div style="font-size:17px;line-height:1.4;">'
               f"<b>{row['ROCKNAME']}</b><br>"
               f"~{int(row['nominal age'])} Ma · Grade {grade}</div>")
    # Only published notebooks (stem prefix <= PUBLISH_AGE_MAX) get a link; older
    # poles are still shown but without one.
    notebook_line = (
        f"<a href='{href}' target='_blank' rel='noopener'>Open notebook &rarr;</a>"
        if row["_publish"] else "")
    # blank for pre-assembly poles from blocks other than the Superior craton
    paleolat_line = (
        "" if pd.isna(row["_duluth_paleolat"])
        else f"Duluth paleolat: {fmt(row['_duluth_paleolat'])}&deg;<br>")
    # sedimentary units are shown at their inclination-shallowing-corrected
    # (Kent) position; give the as-measured pole beneath it
    if row.get("_kent", False):
        pole_lines = (
            f"Pole (incl.-corrected): {fmt(row['PLAT'])}&deg;N, "
            f"{fmt(row['PLONG'])}&deg;E (A95 {fmt(row['A95'])}&deg;)<br>"
            f'<span style="color:#666;">as measured: '
            f"{fmt(row['_plat_uncorrected'])}&deg;N, "
            f"{fmt(row['_plong_uncorrected'])}&deg;E "
            f"(A95 {fmt(row['_a95_uncorrected'])}&deg;)</span><br>")
    else:
        pole_lines = (f"Pole: {fmt(row['PLAT'])}&deg;N, "
                      f"{fmt(row['PLONG'])}&deg;E "
                      f"(A95 {fmt(row['A95'])}&deg;)<br>")
    popup_html = (
        f'<div style="font-size:16px;line-height:1.5;">'
        f"<b>{row['ROCKNAME']}</b><br>"
        f"<i>{row['Terrane']}</i><br>"
        f"~{int(row['nominal age'])} Ma &middot; Grade {grade}<br>"
        f"Site: {fmt(row['SLAT'])}&deg;N, {fmt(slon)}&deg;E<br>"
        f"{pole_lines}"
        f"{paleolat_line}"
        f"{short_reference(row['POLE AUTHORS'], row['YEAR'])}"
        + (f"<br>{notebook_line}" if notebook_line else "")
        + "</div>")
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
    # Base map, but control=False keeps it out of the layer toggle.
    folium.TileLayer("CartoDB positron", control=False).add_to(m)

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
        # Published notebooks (stem prefix <= PUBLISH_AGE_MAX) link; older poles
        # are still listed as plain text.
        link = (f"[{row['ROCKNAME']}](pole_notebooks/{row['_stem']}.ipynb)"
                if row["_publish"] else str(row["ROCKNAME"]))
        # Every pole is listed as measured. A dagger marks the sedimentary units
        # whose inclination-shallowing correction is tabulated separately below.
        if row.get("_kent", False):
            link += "&dagger;"
        lines.append(
            f"| {row['Terrane']} | {link} | {int(row['nominal age'])} "
            f"| {grade_label(row['Grade'])} "
            f"| {fmt(row['SLONG'] % 360.0)} | {fmt(row['SLAT'])} "
            f"| {fmt(row['_plong_uncorrected'])} "
            f"| {fmt(row['_plat_uncorrected'])} "
            f"| {fmt(row['_a95_uncorrected'])} "
            f"| {fmt(row['_duluth_paleolat_uncorrected'])} "
            f"| {short_reference(row['POLE AUTHORS'], row['YEAR'])} |"
        )
    return "\n".join(lines)


def build_kent_table(df):
    """Build the Markdown inclination-shallowing table for the landing page.

    A simplified companion to ``data/nordic_summaries/kent_pole_table.tex``:
    the flattening factor, the corrected (Kent) mean pole and the semi-angles of
    its 95% confidence ellipse, and the Duluth paleolatitude before and after
    correction. The ellipse axis *directions* carried in the LaTeX table are
    dropped here, since the axes are not independently useful without them being
    plotted -- ``kent_poles_combined.csv`` has the full parameterization.

    Args:
        df (pd.DataFrame): Output of :func:`load_poles`.

    Returns:
        str: A GitHub-flavored Markdown table, sedimentary units only.
    """
    kent = load_kent_poles()
    header = ("| Unit | Age (Ma) | *f* source | *f* | Plon | Plat "
              "| &zeta;95 | &eta;95 | Duluth paleolat as measured "
              "| Duluth paleolat corrected |")
    sep = "|" + "|".join(["---"] * 10) + "|"
    lines = [header, sep]
    for _, row in df[df["_kent"]].iterrows():
        src = kent.get(str(row["ROCKNAME"]).strip(), {})
        f_val, f_lo, f_hi = (src.get("f"), src.get("f low"), src.get("f high"))
        f_txt = (f"{float(f_val):.2f} ({float(f_lo):.2f}&ndash;{float(f_hi):.2f})"
                 if f_val not in (None, "") else "")
        link = (f"[{row['ROCKNAME']}](pole_notebooks/{row['_stem']}.ipynb)"
                if row["_publish"] else str(row["ROCKNAME"]))
        lines.append(
            f"| {link} | {int(row['nominal age'])} "
            f"| {row['_kent_method']} | {f_txt} "
            f"| {fmt(row['PLONG'])} | {fmt(row['PLAT'])} "
            f"| {fmt(row['_Zeta'])} | {fmt(row['_Eta'])} "
            f"| {fmt(row['_duluth_paleolat_uncorrected'])} "
            f"| {fmt(row['_duluth_paleolat'])} |"
        )
    return "\n".join(lines)


USE_COLS = ["Terrane", "ROCKNAME", "SLAT", "SLONG", "PLAT", "PLONG", "A95",
            "nominal age", "Grade", "POLE AUTHORS", "YEAR"]


def _prep(df):
    """Coerce ``nominal age`` to int and drop rows lacking a plottable position."""
    df = df.copy()
    df["nominal age"] = pd.to_numeric(df["nominal age"], errors="coerce")
    df = df.dropna(subset=["SLAT", "SLONG", "PLAT", "PLONG", "nominal age"])
    df["nominal age"] = df["nominal age"].astype(int)
    return df


def load_poles():
    """Assemble the table/map rows: site-level recreated poles + older poles.

    Poles at or below ``PUBLISH_AGE_MAX`` are taken from the per-notebook
    Nordic summaries (``COMBINED_CSV``), recreated at the site level and each
    linked to its assessment notebook. This is the set of poles the
    accompanying manuscript tabulates, and the same source that
    ``data/nordic_summaries/combine_nordic_summaries.py`` writes the manuscript
    table from, so the two agree row for row over that interval.

    Older poles (> ``PUBLISH_AGE_MAX``) are carried over from previous
    compilations via ``CSV``. They have not been recreated at the site level
    and are listed without a notebook link, but they are retained here because
    the compilation covers Laurentia and its constituent blocks back to the
    Archean even though the manuscript does not discuss them. The two age
    ranges are disjoint at ``PUBLISH_AGE_MAX``, so the union introduces no
    duplicates.

    ``_duluth_paleolat`` is computed from the pole after rotating rifted
    terranes into Laurentia coordinates via :data:`TERRANE_EULER_POLES`, and is
    NaN for poles older than ``PUBLISH_AGE_MAX`` whose terrane is not in
    :data:`DULUTH_PALEOLAT_TERRANES`.

    Returns:
        pd.DataFrame: rows with ``_stem`` (notebook stem), ``_publish`` (whether
        a notebook link is emitted) and ``_duluth_paleolat``.
    """
    # site-level recreated poles (<= PUBLISH_AGE_MAX) from the Nordic summaries
    summ = _prep(pd.read_csv(COMBINED_CSV))
    stem_by_rock = summary_stem_by_rockname()
    summ["_stem"] = summ["ROCKNAME"].map(stem_by_rock)
    unresolved = summ["_stem"].isna()
    if unresolved.any():
        summ.loc[unresolved, "_stem"] = [
            convention_stem(rn, ag) for rn, ag in
            zip(summ.loc[unresolved, "ROCKNAME"], summ.loc[unresolved, "nominal age"])]

    # older poles (> PUBLISH_AGE_MAX) carried over from previous compilations
    legacy = _prep(pd.read_csv(CSV))
    legacy = legacy[legacy["nominal age"] > PUBLISH_AGE_MAX].copy()
    # a pole recreated at the site level supersedes its carried-over counterpart,
    # which may sit just the other side of PUBLISH_AGE_MAX with a slightly
    # different nominal age (e.g. the ECMB dykes, recreated at 1779 Ma and
    # carried over at 1780 Ma) and would otherwise be listed twice
    superseded = legacy["ROCKNAME"].isin(set(summ["ROCKNAME"]))
    if superseded.any():
        print("  Legacy rows superseded by a site-level recreation: "
              + ", ".join(sorted(legacy.loc[superseded, "ROCKNAME"])))
        legacy = legacy[~superseded].copy()
    legacy["_stem"] = legacy["ROCKNAME"].map(build_stem_map(legacy))

    df = pd.concat([summ[USE_COLS + ["_stem"]], legacy[USE_COLS + ["_stem"]]],
                   ignore_index=True)
    # The map plots sedimentary units at their inclination-shallowing-corrected
    # (Kent) position; ``apply_kent_poles`` keeps the as-measured values in the
    # ``_*_uncorrected`` columns so the compilation table can list those instead.
    df = apply_kent_poles(df)
    df["_publish"] = df["_stem"].map(is_published)

    def _duluth(lat_col, lon_col):
        """Duluth paleolatitude from a pole, rotating rifted terranes first."""
        rot_lat = df[lat_col].astype(float).copy()
        rot_lon = df[lon_col].astype(float).copy()
        for idx in df.index:
            euler = TERRANE_EULER_POLES.get(str(df.at[idx, "Terrane"]))
            if euler is None or pd.isna(rot_lat[idx]) or pd.isna(rot_lon[idx]):
                continue
            rlat, rlon = pmag.pt_rot(euler, [rot_lat[idx]], [rot_lon[idx]])
            rot_lat[idx], rot_lon[idx] = rlat[0], rlon[0]
        return paleolatitude(DULUTH_LAT, DULUTH_LON, rot_lat, rot_lon)

    df["_duluth_paleolat"] = _duluth("PLAT", "PLONG")
    df["_duluth_paleolat_uncorrected"] = _duluth("_plat_uncorrected",
                                                 "_plong_uncorrected")
    # before Laurentia assembled, only Superior/Trans-Hudson poles constrain Duluth
    pre_assembly = ((df["nominal age"] > PUBLISH_AGE_MAX)
                    & ~df["Terrane"].astype(str).isin(DULUTH_PALEOLAT_TERRANES))
    df.loc[pre_assembly, "_duluth_paleolat"] = np.nan
    df.loc[pre_assembly, "_duluth_paleolat_uncorrected"] = np.nan
    return df.sort_values("nominal age", kind="stable").reset_index(drop=True)


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
and zoom/pan freely.

Poles at or below 1779 Ma are the site-level recreations compiled in \
`data/nordic_summaries/`, each with an assessment notebook; older poles are \
carried over from `data/Laurentia_poles.csv`, the compilation assembled at \
earlier Nordic Paleomagnetic Workshops, and are shown without a link. Older \
than 1780 Ma the Duluth paleolatitude is only given for poles from the \
Superior craton and the Trans-Hudson orogen, since a pole from a block that \
had not yet joined Laurentia does not constrain where Duluth was.

Detrital remanence is shallowed during compaction, so the sedimentary units \
are shown at their inclination-shallowing-corrected (Kent mean) position \
rather than as measured; their popups give the as-measured pole beneath it. \
The corrections are tabulated in \
`data/nordic_summaries/kent_poles_combined.csv`."""

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


def write_table_into_compilation(table_md, kent_table_md=None):
    """Replace the table blocks in compilation.md between the sentinel comments.

    Args:
        table_md (str): The main pole table, written between
            :data:`TABLE_START` / :data:`TABLE_END`.
        kent_table_md (str | None): The inclination-shallowing table, written
            between :data:`KENT_TABLE_START` / :data:`KENT_TABLE_END`. Skipped
            when None.
    """
    with open(COMPILATION_MD, encoding="utf-8") as fh:
        text = fh.read()
    blocks = [(TABLE_START, TABLE_END, table_md)]
    if kent_table_md is not None:
        blocks.append((KENT_TABLE_START, KENT_TABLE_END, kent_table_md))
    for start, end, md in blocks:
        if start not in text or end not in text:
            raise SystemExit(
                f"Sentinels {start}/{end} not found in compilation.md; "
                "add them where the table should go.")
        text = re.sub(re.escape(start) + r".*?" + re.escape(end),
                      f"{start}\n{md}\n{end}", text, flags=re.DOTALL)
    with open(COMPILATION_MD, "w", encoding="utf-8") as fh:
        fh.write(text)


def report_link_coverage(df):
    """Report how pole links map onto the notebooks that actually exist.

    Flags two kinds of drift so naming mismatches surface at build time:
    poles whose linked stem has no file on disk, and notebooks present in
    ``pole_notebooks/`` that no pole row links to.
    """
    nb_dir = os.path.join(ROOT, "pole_notebooks")
    on_disk = {f[:-6] for f in os.listdir(nb_dir) if f.endswith(".ipynb")}
    # only the published (linked) poles need a resolving notebook
    published = df[df["_publish"]]
    linked = {stem: rn for stem, rn in zip(published["_stem"],
                                           published["ROCKNAME"])}

    resolved = sum(stem in on_disk for stem in linked)
    print(f"{resolved}/{len(linked)} published pole links resolve to a notebook.")

    missing = [(rn, stem) for stem, rn in linked.items() if stem not in on_disk]
    if missing:
        print("  WARNING: published links with no notebook on disk:")
        for rn, st in sorted(missing):
            print(f"    {rn!r} -> pole_notebooks/{st}.ipynb")

    # published notebooks (stem prefix <= PUBLISH_AGE_MAX) that no pole links to
    orphans = sorted(st for st in on_disk - set(linked) if is_published(st))
    if orphans:
        print("  Published notebooks on disk that no pole links to (check naming):")
        for st in orphans:
            print(f"    pole_notebooks/{st}.ipynb")


def main():
    df = load_poles()

    os.makedirs(os.path.dirname(MAP_HTML), exist_ok=True)
    build_map(df).save(MAP_HTML)          # standalone copy for direct preview
    write_pole_map_notebook(df)
    write_table_into_compilation(build_table(df), build_kent_table(df))

    print(f"Wrote {os.path.relpath(MAP_HTML, ROOT)} (preview, {len(df)} poles)")
    print(f"Wrote {os.path.relpath(POLE_MAP_IPYNB, ROOT)} (interactive map)")
    print(f"Updated table block in {os.path.relpath(COMPILATION_MD, ROOT)} "
          f"({len(df)} rows)")
    print(f"Map notebook links use BASE_URL={BASE_URL!r}")
    report_link_coverage(df)


if __name__ == "__main__":
    main()
