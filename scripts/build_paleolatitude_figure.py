"""Two-panel Duluth-paleolatitude figures comparing the 2017 and updated poles.

Recreates the top panel of Figure 4.5 of Swanson-Hysell (2021, The Precambrian
paleogeography of Laurentia) as two stacked panels:

- Upper: the compilation as of the 2017 Nordic workshop
  (``data/Evans_et_al_2021_compilation.csv``).
- Lower: the updated compilation from the 2022 (Kringerdalen) and 2026
  (Iloranta) workshops, taken from the recreated-from-site-level summaries
  (``data/nordic_summaries/nordic_summaries_combined.csv``).

Each marker is the paleolatitude of Duluth, Minnesota (lat 46.79N, lon 92.10W)
implied by a pole, with vertical error bars from the pole's confidence region
and horizontal error bars from the magmin/magmax age bounds, colored/marked by
terrane and reliability grade. Poles whose age half-range reaches 50 Myr are
excluded, as is Svalbard (no longer considered Laurentia sensu stricto). A
geologic-timescale strip (eras over periods) anchors the base.

In the updated panel the sedimentary units are plotted at their
inclination-shallowing-corrected (Kent mean) positions from
``data/nordic_summaries/kent_poles_combined.csv``, and their error bar is the
Kent ellipse projected onto the pole-to-Duluth great circle (see
:func:`paleolat_uncertainty`) rather than its equal-area circular
approximation, since the flattening-factor uncertainty elongates the ellipse
along that direction; the 2017 panel is left as Evans et al. (2021) published
it, so the two panels also show what correcting for compaction shallowing does
to the record.

Two versions are written to ``_static/``: a full 1800-700 Ma view and a
1180-700 Ma view (the Keweenawan-Grenville interval)::

    python scripts/build_paleolatitude_figure.py
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pmagpy.pmag as pmag  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_pole_map import (  # noqa: E402
    DULUTH_LAT, DULUTH_LON, KENT_ELLIPSE_COLUMNS, ROOT, apply_kent_poles,
    paleolatitude)

# Myr; poles whose age *half-range* ((himagage - lomagage) / 2) reaches this are
# excluded. The test is strict (half-range < MAX_AGE_UNCERTAINTY) and is on the
# half-range rather than the one-sided distances from the nominal age, which
# matters for asymmetric ranges: Torridon (975 Ma, 950-1050) is +75/-25 with a
# half-range of exactly 50, so it fails the vetting and enters these figures only
# through include_scotland. Note that build_apwp_figure.py vets the SphereUDE
# path fit with a non-strict test (half_unc > FIT_MAX_AGE_HALF_UNC drops), so
# Torridon is retained in the path fit while being excluded here.
MAX_AGE_UNCERTAINTY = 50
EVANS_CSV = os.path.join(ROOT, "data", "Evans_et_al_2021_compilation.csv")
# Updated compilation = the recreated-from-site-level poles emitted by the
# notebooks and merged by combine_nordic_summaries.py (grows as notebooks are
# built), including the rebuilt Greenland poles.
UPDATED_CSV = os.path.join(ROOT, "data", "nordic_summaries",
                           "nordic_summaries_combined.csv")
OUT_BASE = os.path.join(ROOT, "_static", "Laurentia_paleolatitude_comparison")

# Euler poles [pole_lat, pole_lon, rotation angle] used to rotate poles from
# terranes separated from Laurentia back into the Laurentia reference frame
# before the Duluth paleolatitude is computed; these mirror
# pole_tools.TERRANE_EULER_POLES (Greenland: Roest & Srivastava, 1989; Scotland:
# Torsvik & Cocks, 2017). The rebuilt Greenland and Scotland poles carry the
# correct 'Laurentia-Greenland' / 'Laurentia-Scotland' terrane labels, in
# present-day coordinates, so they are rotated here. Svalbard is excluded from
# the figure.
TERRANE_EULER_POLES = {
    "Laurentia-Greenland": [67.5, -118.5, -13.8],
    "Laurentia-Greenland-Nain": [67.5, -118.5, -13.8],
    "Laurentia-Scotland": [78.6, 161.9, -31.0],
}

# (age_min, age_max) for each version produced.
VERSIONS = [(700, 1800), (700, 1180)]

# (label, color, marker), in legend order. Markers/colors follow the Fig 4.5
# scheme: Laurentia split by grade, associated terranes by region.
CATEGORIES = [
    ("Laurentia A poles", "#034C8C", "D"),
    ("Laurentia B poles", "#0477BF", "s"),
    ("Trans-Hudson orogen poles", "#4A90D9", "^"),
    ("Greenland poles", "#808000", "p"),
    ("Scotland poles", "#C8A415", "v"),
    ("Other Laurentia terrane", "#888888", "d"),
]

# Poles that fail the age vetting are drawn as open (hollow) markers to set them
# apart from the age-vetted poles. In practice these are the two Scotland poles,
# Stoer (+/-70 Myr) and Torridon (half-range exactly 50 Myr), retained via
# include_scotland in the "with Scotland" variants only.
# Terrane is still encoded by marker shape and color, independent of vetting.
SCOTLAND_LABEL = "Scotland poles"


def marker_style(vetted, color):
    """Per-pole face/edge styling for the errorbar markers.

    Poles that passed the age vetting are filled; those retained despite failing
    it (via ``include_scotland``) are hollow.
    """
    if not vetted:
        return dict(markerfacecolor="none", markeredgecolor=color,
                    markeredgewidth=1.2)
    return dict(markerfacecolor=color, markeredgecolor="black",
                markeredgewidth=0.4)

# Proterozoic eras and periods: (name, top_age, base_age) in Ma.
TIMESCALE_ERAS = [
    ("Paleoproterozoic", 1600, 2500),
    ("Mesoproterozoic", 1000, 1600),
    ("Neoproterozoic", 538.8, 1000),
]
TIMESCALE_PERIODS = [
    ("Statherian", 1600, 1800),
    ("Calymmian", 1400, 1600),
    ("Ectasian", 1200, 1400),
    ("Stenian", 1000, 1200),
    ("Tonian", 720, 1000),
    ("Cryogenian", 635, 720),
]


# Full-height tectonic-phase bands of the Rodinia cycle:
# (label, young_age, old_age, color); young_age None means the panel's young
# edge.
PHASES = [
    ("Grenvillian Orogeny\n(Rodinia assembly)", 985, 1085, "#d1564f"),
    ("stable supercontinent\ninterior (Rodinia)", 780, 985, "#4e9a51"),
    ("rifting\n(Rodinia\nbreak-up)", None, 780, "#dbb23a"),
]


def categorize(terrane, grade):
    """Map (Terrane, Grade) to a (label, color, marker) plotting category."""
    terrane = str(terrane)
    if terrane == "Laurentia":
        return CATEGORIES[0] if grade == "A" else CATEGORIES[1]
    if "Trans-Hudson" in terrane:
        return CATEGORIES[2]
    if "Greenland" in terrane:
        return CATEGORIES[3]
    if "Scotland" in terrane:
        return CATEGORIES[4]
    return CATEGORIES[5]


def rotate_to_laurentia(d):
    """Rotate separated-terrane (Greenland, Scotland) poles into the Laurentia
    reference frame in place, so the Duluth paleolatitude is computed in
    Laurentia coordinates. Poles of any terrane without an entry in
    ``TERRANE_EULER_POLES`` are left unchanged.

    Where a Kent ellipse is carried alongside the pole (the ``_Zdec``/``_Zinc``
    and ``_Edec``/``_Einc`` axis directions added by ``apply_kent_poles``) its
    axes are rotated by the same Euler rotation. A rigid rotation carries the
    ellipse with the pole, so the axis directions have to move with it for
    :func:`paleolat_uncertainty` to project the ellipse correctly.
    """
    axis_pairs = [("_Zinc", "_Zdec"), ("_Einc", "_Edec")]
    for idx in d.index:
        euler = TERRANE_EULER_POLES.get(str(d.at[idx, "Terrane"]))
        if euler is None:
            continue
        plat, plon = d.at[idx, "PLAT"], d.at[idx, "PLONG"]
        if pd.isna(plat) or pd.isna(plon):
            continue
        rlat, rlon = pmag.pt_rot(euler, [plat], [plon])
        d.at[idx, "PLAT"], d.at[idx, "PLONG"] = rlat[0], rlon[0]
        for lat_col, lon_col in axis_pairs:
            if lat_col not in d.columns or pd.isna(d.at[idx, lat_col]):
                continue
            alat, alon = pmag.pt_rot(euler, [d.at[idx, lat_col]],
                                     [d.at[idx, lon_col]])
            d.at[idx, lat_col], d.at[idx, lon_col] = alat[0], alon[0]
    return d


def bearing(lat1, lon1, lat2, lon2):
    """Initial great-circle bearing from point 1 to point 2, degrees E of N."""
    phi1, phi2 = np.radians([float(lat1), float(lat2)])
    dlon = np.radians(float(lon2) - float(lon1))
    return np.degrees(np.arctan2(
        np.sin(dlon) * np.cos(phi2),
        np.cos(phi1) * np.sin(phi2)
        - np.sin(phi1) * np.cos(phi2) * np.cos(dlon))) % 360.0


def paleolat_uncertainty(d, site_lat=DULUTH_LAT, site_lon=DULUTH_LON):
    """Uncertainty in the paleolatitude each pole implies for the site, in degrees.

    The paleolatitude is 90 minus the angular distance from the site to the
    pole, so an uncertainty in pole position maps one-to-one onto paleolatitude
    along the site-to-pole great circle. For a pole with a circular A95 that is
    just A95, whatever the pole's orientation. For the inclination-shallowing-
    corrected (Kent) poles of the sedimentary units the confidence region is an
    ellipse elongated along the direction the shallowing correction displaces
    the pole, which is close to the site-to-pole great circle, so the relevant
    radius is the ellipse evaluated in that direction,

        r(theta) = ((cos(theta) / zeta95)^2 + (sin(theta) / eta95)^2)^(-1/2),

    with ``theta`` the angle between the ellipse's major axis and the bearing
    from the pole to the site. Using the equal-area circular radius
    ``sqrt(zeta95 * eta95)`` that ``A95`` carries for these rows instead would
    understate the paleolatitude uncertainty, because zeta95 is the semi-axis
    that is nearly aligned with the site.

    Args:
        d (pd.DataFrame): Poles with ``PLAT``/``PLONG``/``A95``, optionally
            carrying the Kent ellipse columns added by ``apply_kent_poles``
            and already rotated into the Laurentia frame.
        site_lat (float): Site latitude (degrees N).
        site_lon (float): Site longitude (degrees E).

    Returns:
        pd.Series: The paleolatitude uncertainty for each row, falling back to
        ``A95`` wherever no ellipse is carried.
    """
    err = pd.to_numeric(d["A95"], errors="coerce")
    if not all("_" + c in d.columns for c in KENT_ELLIPSE_COLUMNS):
        return err
    for idx in d.index:
        zeta, eta = d.at[idx, "_Zeta"], d.at[idx, "_Eta"]
        if pd.isna(zeta) or pd.isna(eta) or zeta <= 0 or eta <= 0:
            continue
        to_site = bearing(d.at[idx, "PLAT"], d.at[idx, "PLONG"],
                          site_lat, site_lon)
        to_major = bearing(d.at[idx, "PLAT"], d.at[idx, "PLONG"],
                           d.at[idx, "_Zinc"], d.at[idx, "_Zdec"])
        theta = np.radians(to_site - to_major)
        err.at[idx] = 1.0 / np.hypot(np.cos(theta) / zeta,
                                     np.sin(theta) / eta)
    return err


def load(path, age_min, age_max, include_scotland=False, kent=False):
    """Load a compilation CSV, restricted to Laurentia poles in the interval.

    By default poles whose age half-range reaches ``MAX_AGE_UNCERTAINTY`` Myr are
    excluded, which drops both Scotland poles: Stoer (+/-70 Myr) and Torridon
    (half-range exactly 50 Myr). If ``include_scotland`` is True, Scotland poles
    are retained regardless of age control, which adds the two of them back as
    lower-confidence points.

    With ``kent=True`` the sedimentary units are replaced by their
    inclination-shallowing-corrected (Kent mean) poles via
    :func:`build_pole_map.apply_kent_poles`, so the paleolatitude plotted for
    them is the corrected one rather than the compaction-shallowed minimum.
    This is applied to the updated compilation only; the 2017 panel is drawn
    from the poles as Evans et al. (2021) published them.
    """
    d = pd.read_csv(path)
    d["Grade"] = d["Grade"].astype(str).str.strip()
    for c in ["PLAT", "PLONG", "A95", "nominal age", "lomagage", "himagage"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    if kent:
        d = apply_kent_poles(d)
    d = d[d["Terrane"].astype(str).str.startswith("Laurentia")]
    # Svalbard is no longer considered Laurentia sensu stricto.
    d = d[~d["Terrane"].astype(str).str.contains("Svalbard")]
    d = d[(d["nominal age"] >= age_min) & (d["nominal age"] <= age_max)].copy()

    # Rotate Greenland/Scotland poles into the Laurentia frame before computing
    # the Duluth paleolatitude.
    d = rotate_to_laurentia(d)
    d["Duluth_plat"] = paleolatitude(DULUTH_LAT, DULUTH_LON, d["PLAT"],
                                     d["PLONG"])
    d["age_hi"] = (d["himagage"] - d["nominal age"]).clip(lower=0).fillna(0)
    d["age_lo"] = (d["nominal age"] - d["lomagage"]).clip(lower=0).fillna(0)
    d["A95"] = d["A95"].fillna(0)
    # Vertical error bar: A95 for the circular poles, and for the Kent poles the
    # ellipse projected onto the pole-to-Duluth great circle rather than its
    # equal-area circular approximation.
    d["plat_err"] = paleolat_uncertainty(d).fillna(0)

    # Drop poles whose age half-range reaches MAX_AGE_UNCERTAINTY Myr, optionally
    # exempting Scotland poles so they can be shown as lower-confidence additions.
    # Poles with a missing lomagage/himagage fall through as 0 and are retained,
    # as before.
    half_unc = ((d["himagage"] - d["lomagage"]) / 2).clip(lower=0).fillna(0)
    vetted = half_unc < MAX_AGE_UNCERTAINTY
    within = vetted
    if include_scotland:
        within = within | d["Terrane"].astype(str).str.contains("Scotland")
    d = d[within].copy()
    # Records which retained poles actually passed the age vetting; those that
    # did not (kept only via include_scotland) are drawn hollow by marker_style.
    d["age_vetted"] = vetted.loc[d.index]

    cats = [categorize(t, g) for t, g in zip(d["Terrane"], d["Grade"])]
    d["cat_label"] = [c[0] for c in cats]
    return d


def load_updated(age_min, age_max, include_scotland=False):
    """Updated-panel poles: the recreated-from-site-level summaries, including
    the rebuilt Greenland poles (now in the combined summaries), each rotated
    into the Laurentia reference frame by :func:`load`.

    The sedimentary units are taken at their inclination-shallowing-corrected
    (Kent mean) position, so the paleolatitudes plotted for them are corrected
    rather than the compaction-shallowed minima. The 2017 panel is left as
    Evans et al. (2021) published it."""
    return load(UPDATED_CSV, age_min, age_max,
                include_scotland=include_scotland, kent=True)


def plot_panel(ax, d, panel_label, age_min, age_max, grenville_xy=(1010, -52)):
    """Plot one Duluth-paleolatitude-versus-age panel."""
    # Full-height transparent tectonic-phase bands, bold-labeled at the top.
    # A wide (compressed) view crowds the bands, so use a smaller font and drop
    # the middle label to a second tier.
    stagger = (age_max - age_min) > 700
    fs = 9 if stagger else 11
    for i, (name, young, old, fill) in enumerate(PHASES):
        lo = max(age_min if young is None else young, age_min)
        hi = min(old, age_max)
        if hi <= lo:
            continue
        ax.axvspan(lo, hi, color=fill, alpha=0.14, lw=0, zorder=0)
        y = 54 if (stagger and i == 1) else 73
        ax.text((lo + hi) / 2, y, name, ha="center", va="top", ma="center",
                fontsize=fs, fontweight="bold", color="#333")

    for label, color, marker in CATEGORIES:
        sub = d[d["cat_label"] == label]
        if sub.empty:
            continue
        # Split the category so vetted poles draw filled and unvetted ones
        # hollow; only the first non-empty part carries the legend entry so the
        # category still appears exactly once.
        labeled = False
        for vetted in (True, False):
            part = sub[sub["age_vetted"] == vetted]
            if part.empty:
                continue
            ax.errorbar(
                part["nominal age"], part["Duluth_plat"],
                yerr=part["plat_err"], xerr=[part["age_lo"], part["age_hi"]],
                fmt=marker, color=color, markersize=6, elinewidth=1,
                capsize=2, linestyle="none",
                label=None if labeled else label,
                **marker_style(vetted, color))
            labeled = True

    ax.axhline(0, color="grey", lw=0.8, linestyle=":", zorder=0)
    if age_min < 1600 < age_max:  # Paleo|Meso boundary only
        ax.axvline(1600, color="lightgrey", lw=0.8, zorder=0)
    ax.text(1135, 60, "Logan\nLoop", fontsize=10, ha="left", style="italic",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.6))
    ax.text(1075, 20, "Keweenawan\nTrack", fontsize=10, style="italic")
    ax.text(grenville_xy[0], grenville_xy[1], "Grenville\nLoop", fontsize=10,
            ha="right", style="italic")

    ax.set_xlim(age_max, age_min)  # older on the left, younger on the right
    ax.set_ylim(-75, 75)
    ax.set_ylabel("Paleolatitude (°)")
    ax.text(0.5, 0.04, panel_label, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="bottom", ha="center",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.7))


def draw_timescale(ax, age_min, age_max):
    """Draw an era (upper row) and period (lower row) timescale strip."""
    def boxes(items, y0, height, fontsize):
        for name, top, base in items:
            x0, x1 = max(top, age_min), min(base, age_max)
            if x1 <= x0:
                continue
            ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, height,
                                       facecolor="white", edgecolor="black",
                                       lw=0.6))
            if x1 - x0 >= 0.06 * (age_max - age_min):  # skip labels on slivers
                ax.text((x0 + x1) / 2, y0 + height / 2, name, ha="center",
                        va="center", fontsize=fontsize)

    boxes(TIMESCALE_ERAS, 0.5, 0.5, 8.5)
    boxes(TIMESCALE_PERIODS, 0.0, 0.5, 7.5)
    ax.set_xlim(age_max, age_min)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("Age (Ma)")


def make_figure(age_min, age_max, out_path, include_scotland=False):
    """Build and save one two-panel comparison figure for an age interval."""
    evans = load(EVANS_CSV, age_min, age_max, include_scotland=include_scotland)
    updated = load_updated(age_min, age_max, include_scotland=include_scotland)

    fig, (ax0, ax1, axts) = plt.subplots(
        3, 1, figsize=(11, 9.3), sharex=True,
        gridspec_kw={"height_ratios": [1, 1, 0.16], "hspace": 0.08})
    # Grenville Loop labels: same paleolatitude in both panels, each just to
    # the left (older side) of that panel's Grenville-age pole (Haliburton
    # ~1015 Ma in 2017; Adirondack ~887 Ma in the update).
    plot_panel(ax0, evans, "2017 compilation (Evans et al., 2021)",
               age_min, age_max, grenville_xy=(1040, -44))
    plot_panel(ax1, updated, "Updated compilation (2022 & 2026 workshops)",
               age_min, age_max, grenville_xy=(915, -44))
    ax0.tick_params(labelbottom=False)
    # Drop the (unlabeled) bottom ticks on the lower data panel so the timescale
    # strip can sit flush beneath it without a redundant tick row.
    ax1.tick_params(labelbottom=False, bottom=False)
    draw_timescale(axts, age_min, age_max)

    # Slide the timescale strip up so its top edge is flush with the bottom of
    # the lower data panel (the gridspec hspace leaves a gap otherwise).
    pos_lower = ax1.get_position()
    pos_ts = axts.get_position()
    axts.set_position([pos_ts.x0, pos_lower.y0 - pos_ts.height,
                       pos_ts.width, pos_ts.height])

    # Place the title just above the top panel rather than floating near the
    # figure edge.
    fig.suptitle("Position of Duluth (lat = 46.79°N, lon = 92.10°W) implied by "
                 "Laurentia poles", fontsize=13,
                 y=ax0.get_position().y1 + 0.02)

    # single legend from the union of categories present in either panel
    handles = {}
    for ax in (ax0, ax1):
        for h, lab in zip(*ax.get_legend_handles_labels()):
            handles.setdefault(lab, h)
    order = [lab for lab, _, _ in CATEGORIES if lab in handles]
    ax0.legend([handles[lab] for lab in order], order, loc="lower left",
               fontsize=9, framealpha=0.9)

    fig.savefig(out_path + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(out_path + ".pdf", bbox_inches="tight")
    plt.close(fig)
    return len(evans), len(updated)


def main():
    os.makedirs(os.path.dirname(OUT_BASE), exist_ok=True)
    for age_min, age_max in VERSIONS:
        # Each interval is produced twice: the default (age-vetted) version and
        # a "_withScotland" version that also shows the loose-age Scotland poles.
        for include_scotland in (False, True):
            suffix = "_withScotland" if include_scotland else ""
            out = f"{OUT_BASE}_{age_max}_{age_min}{suffix}"
            n_evans, n_updated = make_figure(age_min, age_max, out,
                                             include_scotland=include_scotland)
            print(f"Wrote {os.path.relpath(out, ROOT)}.png / .pdf "
                  f"({age_max}-{age_min} Ma: Evans {n_evans}, "
                  f"updated {n_updated} poles)")


if __name__ == "__main__":
    main()
