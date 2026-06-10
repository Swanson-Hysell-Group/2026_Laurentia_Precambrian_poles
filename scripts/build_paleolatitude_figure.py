"""Two-panel Duluth-paleolatitude figures comparing the 2017 and updated poles.

Recreates the top panel of Figure 4.5 of Swanson-Hysell (2021, The Precambrian
paleogeography of Laurentia) as two stacked panels:

- Upper: the compilation as of the 2017 Nordic workshop
  (``data/Evans_et_al_2021_compilation.csv``).
- Lower: the updated compilation from the 2022 (Kringerdalen) and 2026
  (Iloranta) workshops, taken from the recreated-from-site-level summaries
  (``data/nordic_summaries/nordic_summaries_combined.csv``).

Each marker is the paleolatitude of Duluth, Minnesota (lat 46.79N, lon 92.10W)
implied by a pole, with vertical error bars from A95 and horizontal error bars
from the magmin/magmax age bounds, colored/marked by terrane and reliability
grade. Poles whose age is uncertain by more than 50 Myr on either bound are
excluded, as is Svalbard (no longer considered Laurentia sensu stricto). A
geologic-timescale strip (eras over periods) anchors the base.

Two versions are written to ``_static/``: a full 1800-700 Ma view and a
1180-700 Ma view (the Keweenawan-Grenville interval)::

    python scripts/build_paleolatitude_figure.py
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_pole_map import (  # noqa: E402
    DULUTH_LAT, DULUTH_LON, ROOT, paleolatitude)

MAX_AGE_UNCERTAINTY = 50  # Myr; poles with a larger +/- age bound are excluded
EVANS_CSV = os.path.join(ROOT, "data", "Evans_et_al_2021_compilation.csv")
# Updated compilation = the recreated-from-site-level poles emitted by the
# notebooks and merged by combine_nordic_summaries.py (grows as notebooks are
# built).
UPDATED_CSV = os.path.join(ROOT, "data", "nordic_summaries",
                           "nordic_summaries_combined.csv")
# Previous compilation, used for now to supply Greenland poles not yet rebuilt.
PREV_CSV = os.path.join(ROOT, "data", "Laurentia_poles.csv")
OUT_BASE = os.path.join(ROOT, "_static", "Laurentia_paleolatitude_comparison")

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


def load(path, age_min, age_max):
    """Load a compilation CSV, restricted to Laurentia poles in the interval."""
    d = pd.read_csv(path)
    d["Grade"] = d["Grade"].astype(str).str.strip()
    for c in ["PLAT", "PLONG", "A95", "nominal age", "lomagage", "himagage"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[d["Terrane"].astype(str).str.startswith("Laurentia")]
    # Svalbard is no longer considered Laurentia sensu stricto.
    d = d[~d["Terrane"].astype(str).str.contains("Svalbard")]
    d = d[(d["nominal age"] >= age_min) & (d["nominal age"] <= age_max)].copy()

    d["Duluth_plat"] = paleolatitude(DULUTH_LAT, DULUTH_LON, d["PLAT"],
                                     d["PLONG"])
    d["age_hi"] = (d["himagage"] - d["nominal age"]).clip(lower=0).fillna(0)
    d["age_lo"] = (d["nominal age"] - d["lomagage"]).clip(lower=0).fillna(0)
    d["A95"] = d["A95"].fillna(0)

    # Drop poles whose age is uncertain by more than MAX_AGE_UNCERTAINTY Myr.
    d = d[(d["age_hi"] <= MAX_AGE_UNCERTAINTY)
          & (d["age_lo"] <= MAX_AGE_UNCERTAINTY)].copy()

    cats = [categorize(t, g) for t, g in zip(d["Terrane"], d["Grade"])]
    d["cat_label"] = [c[0] for c in cats]
    return d


def load_updated(age_min, age_max):
    """Updated-panel poles: the rebuilt summaries, but with Greenland poles
    taken from the previous compilation for now (not all are rebuilt yet)."""
    rebuilt = load(UPDATED_CSV, age_min, age_max)
    rebuilt = rebuilt[~rebuilt["Terrane"].astype(str).str.contains("Greenland")]
    prev = load(PREV_CSV, age_min, age_max)
    greenland = prev[prev["Terrane"].astype(str).str.contains("Greenland")]
    return pd.concat([rebuilt, greenland], ignore_index=True)


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
        ax.errorbar(
            sub["nominal age"], sub["Duluth_plat"],
            yerr=sub["A95"], xerr=[sub["age_lo"], sub["age_hi"]],
            fmt=marker, color=color, markersize=6, elinewidth=1,
            capsize=2, linestyle="none", markeredgecolor="black",
            markeredgewidth=0.4, label=label)

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


def make_figure(age_min, age_max, out_path):
    """Build and save one two-panel comparison figure for an age interval."""
    evans = load(EVANS_CSV, age_min, age_max)
    updated = load_updated(age_min, age_max)

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
    ax1.tick_params(labelbottom=False)
    draw_timescale(axts, age_min, age_max)

    fig.suptitle("Position of Duluth (lat = 46.79°N, lon = 92.10°W) implied by "
                 "Laurentia poles", fontsize=13, y=0.92)

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
        out = f"{OUT_BASE}_{age_max}_{age_min}"
        n_evans, n_updated = make_figure(age_min, age_max, out)
        print(f"Wrote {os.path.relpath(out, ROOT)}.png / .pdf "
              f"({age_max}-{age_min} Ma: Evans {n_evans}, updated {n_updated} "
              f"poles)")


if __name__ == "__main__":
    main()
