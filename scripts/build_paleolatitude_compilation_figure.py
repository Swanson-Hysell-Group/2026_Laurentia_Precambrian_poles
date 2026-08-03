"""Single-panel Duluth-paleolatitude figure for the updated pole compilation.

A companion to ``build_paleolatitude_figure.py`` (which produces the two-panel
2017-versus-updated *comparison*). This script plots only the updated
compilation and extends the record from 1800 Ma to the present day, so the
full Proterozoic-through-Phanerozoic paleolatitude history of Duluth is shown
in a single panel.

Two pole sources are combined:

- Precambrian (700-1800 Ma): the updated, recreated-from-site-level Laurentia
  summaries from the 2022 (Kringerdalen) and 2026 (Iloranta) workshops
  (``data/nordic_summaries/nordic_summaries_combined.csv``), loaded and rotated
  into the Laurentia reference frame by :func:`build_paleolatitude_figure.load`.
- Phanerozoic (0-~540 Ma): the Laurentia poles of the Torsvik et al. (2012)
  global compilation (``data/Torsvik_Laurentia_Pole_Compilation.csv``, the
  Laurentia subset of the GPDB "World701" table). This mirrors the Phanerozoic
  treatment in the Laurentia_Paleogeography ``Laurentia_paleolatitude.pdf``
  figure. These poles are tabulated as south poles, so the GAD paleolatitude
  is negated to place Duluth in its correct (northern, by the late Paleozoic)
  hemisphere.

Each marker is the paleolatitude of Duluth, Minnesota (lat 46.79N, lon 92.10W)
implied by a pole. Full-height bands mark the durations of the three
supercontinents that Laurentia participated in:

- Nuna: 1820-1380 Ma. Assembly is taken as the Trans-Hudson orogeny (ca. 1820
  Ma), which sutured the Archean cratons into the Laurentian core. Breakup is
  taken as ca. 1380 Ma, following Ding et al. (2025): combined with the
  1.45-1.04 Ga paleomagnetic and geologic record from Laurentia, the North
  China Craton, Baltica, and Australia, the divergence of their apparent polar
  wander paths indicates that the core of Nuna (Laurentia, Baltica, Siberia)
  separated from East Nuna (Australia, the NCC) at ca. 1.38 Ga.
- Rodinia: 1060-718 Ma. Assembly is taken as peak Ottawan (Grenvillian)
  metamorphism at ca. 1060 Ma. Breakup is taken as the ca. 718 Ma Franklin
  large igneous province, by which time the conjugates of the western
  Laurentian margin had separated.
- Pangea: 320-175 Ma. Standard assembly-to-onset-of-breakup interval.

Output is written to ``_static/Laurentia_paleolatitude_compilation_1800_0.{png,pdf}``::

    python scripts/build_paleolatitude_compilation_figure.py
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
from build_paleolatitude_figure import (  # noqa: E402
    CATEGORIES, load_updated, marker_style)

TORSVIK_CSV = os.path.join(ROOT, "data",
                           "Torsvik_Laurentia_Pole_Compilation.csv")
OUT_BASE = os.path.join(ROOT, "_static",
                        "Laurentia_paleolatitude_compilation")

# Full age span of the panel: present day to 1800 Ma.
AGE_MIN, AGE_MAX = 0, 1800
# The updated compilation is populated from ca. 1800 Ma forward to ca. 565 Ma
# (Ediacaran); the Phanerozoic from ca. 532 Ma is carried by the Torsvik poles.
# The cutoff sits just above the oldest Torsvik pole so the two records hand off
# without overlap.
PRECAMBRIAN_MIN = 535

# Ediacaran-aged updated poles (younger than the 635 Ma Cryogenian-Ediacaran
# boundary) are drawn at reduced opacity given their larger uncertainties.
EDIACARAN_BASE = 635
EDIACARAN_ALPHA = 0.35

# Torsvik Phanerozoic poles are plotted as their own category.
TORSVIK_LABEL = "Laurentia poles (Torsvik et al., 2012)"
TORSVIK_COLOR = "#8C633F"
TORSVIK_MARKER = "o"

# Two limestone poles excluded in the Laurentia_Paleogeography figure (anomalous
# inclination-shallowed carbonates); excluded here for consistency.
TORSVIK_EXCLUDE = ["St. George  Group  limestone",
                   "Tablehead Group  limestone Mean"]

# Full-height supercontinent-duration bands: (label, young_age, old_age, color).
# See the module docstring for the rationale behind each age bound.
SUPERCONTINENTS = [
    ("Nuna", 1380, 1820, "#4A7BA6"),
    ("Rodinia", 718, 1060, "#4e9a51"),
    ("Pangea", 175, 320, "#d1564f"),
]

# Eras (upper timescale row) and periods (lower row):
# (full_name, abbreviation, top_age, base_age). The longest label that fits
# inside the box is drawn (full name, else abbreviation, else nothing).
TIMESCALE_ERAS = [
    ("Paleoproterozoic", "Paleoprot.", 1600, 2500),
    ("Mesoproterozoic", "Mesoprot.", 1000, 1600),
    ("Neoproterozoic", "Neoprot.", 538.8, 1000),
    ("Paleozoic", "Pz", 251.902, 538.8),
    ("Mesozoic", "Mz", 66.0, 251.902),
    ("Cenozoic", "Cz", 0, 66.0),
]
TIMESCALE_PERIODS = [
    ("Statherian", "Stath.", 1600, 1800),
    ("Calymmian", "Calym.", 1400, 1600),
    ("Ectasian", "Ect.", 1200, 1400),
    ("Stenian", "Sten.", 1000, 1200),
    ("Tonian", "Ton.", 720, 1000),
    ("Cryogenian", "Cryo.", 635, 720),
    ("Ediacaran", "Ediac.", 538.8, 635),
    ("Cambrian", "Cm", 485.4, 538.8),
    ("Ordovician", "O", 443.8, 485.4),
    ("Silurian", "S", 419.2, 443.8),
    ("Devonian", "D", 358.9, 419.2),
    ("Carboniferous", "C", 298.9, 358.9),
    ("Permian", "P", 251.902, 298.9),
    ("Triassic", "Tr", 201.4, 251.902),
    ("Jurassic", "J", 145.0, 201.4),
    ("Cretaceous", "K", 66.0, 145.0),
    ("Paleogene", "Pg", 23.03, 66.0),
    ("Neogene", "Ng", 2.58, 23.03),
    ("Quaternary", "Q", 0, 2.58),
]


def load_torsvik():
    """Phanerozoic Duluth paleolatitudes from the Torsvik et al. (2012)
    Laurentia poles. The compilation tabulates south poles, so the GAD
    paleolatitude is negated to place Duluth in the correct hemisphere."""
    d = pd.read_csv(TORSVIK_CSV)
    d = d[~d["Formation"].isin(TORSVIK_EXCLUDE)].copy()
    for c in ["CLat", "CLon", "A95", "Age"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["Duluth_plat"] = -paleolatitude(DULUTH_LAT, DULUTH_LON, d["CLat"],
                                      d["CLon"])
    return d


def plot_panel(ax, include_scotland=False):
    """Plot the combined Precambrian + Phanerozoic Duluth-paleolatitude panel."""
    updated = load_updated(PRECAMBRIAN_MIN, AGE_MAX,
                           include_scotland=include_scotland)
    torsvik = load_torsvik()

    # Full-height supercontinent-duration bands, bold-labeled at the top.
    for name, young, old, fill in SUPERCONTINENTS:
        lo, hi = max(young, AGE_MIN), min(old, AGE_MAX)
        if hi <= lo:
            continue
        ax.axvspan(lo, hi, color=fill, alpha=0.14, lw=0, zorder=0)
        ax.text((lo + hi) / 2, 73, name, ha="center", va="top",
                fontsize=12, fontweight="bold", color="#333")

    # Updated compilation, split by terrane/grade category. Ediacaran-aged
    # poles (younger than ca. 635 Ma) are drawn at reduced opacity to flag the
    # larger uncertainties associated with the terminal-Neoproterozoic record.
    for label, color, marker in CATEGORIES:
        sub = updated[updated["cat_label"] == label]
        if sub.empty:
            continue
        # Vetted poles draw filled, unvetted ones hollow; only the first
        # non-empty group carries the legend entry so each category appears once.
        labeled = False
        for vetted in (True, False):
            part = sub[sub["age_vetted"] == vetted]
            if part.empty:
                continue
            older = part[part["nominal age"] >= EDIACARAN_BASE]
            ediac = part[part["nominal age"] < EDIACARAN_BASE]
            style = marker_style(vetted, color)
            ax.errorbar(
                older["nominal age"], older["Duluth_plat"],
                yerr=older["plat_err"], xerr=[older["age_lo"], older["age_hi"]],
                fmt=marker, color=color, markersize=6, elinewidth=1,
                capsize=2, linestyle="none",
                label="_nolegend_" if labeled else label, **style)
            labeled = True
            if not ediac.empty:
                ax.errorbar(
                    ediac["nominal age"], ediac["Duluth_plat"],
                    yerr=ediac["plat_err"],
                    xerr=[ediac["age_lo"], ediac["age_hi"]],
                    fmt=marker, color=color, markersize=6, elinewidth=1,
                    capsize=2, linestyle="none", alpha=EDIACARAN_ALPHA,
                    label="_nolegend_", **style)

    # Phanerozoic Torsvik poles.
    ax.errorbar(
        torsvik["Age"], torsvik["Duluth_plat"], yerr=torsvik["A95"],
        fmt=TORSVIK_MARKER, color=TORSVIK_COLOR, markersize=6, elinewidth=1,
        capsize=2, linestyle="none", markeredgecolor="black",
        markeredgewidth=0.4, label=TORSVIK_LABEL, zorder=1000)

    ax.axhline(0, color="grey", lw=0.8, linestyle=":", zorder=0)
    ax.text(1110, 70, "Logan Loop", fontsize=10, ha="center", va="center",
            style="italic",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.6))
    ax.text(1075, 30, "Keweenawan\nTrack", fontsize=10, style="italic")
    ax.text(915, -49, "Grenville\nLoop", fontsize=10, ha="right",
            style="italic")

    ax.set_xlim(AGE_MAX, AGE_MIN)  # older on the left, younger on the right
    ax.set_ylim(-75, 75)
    ax.set_ylabel("Paleolatitude (°)")


def _fits(ax, renderer, text, x_center, x0, x1):
    """True if a rendered text artist is narrower than its box (with margin)."""
    box_px = abs(ax.transData.transform((x1, 0))[0]
                 - ax.transData.transform((x0, 0))[0])
    return text.get_window_extent(renderer=renderer).width <= 0.9 * box_px


def draw_timescale(ax, renderer):
    """Draw an era (upper row) and period (lower row) timescale strip. For each
    box the longest label that fits is drawn: full name, else abbreviation, else
    nothing."""
    # Limits must be set before measuring text, since the fit test relies on
    # transData (which depends on the x-limits).
    ax.set_xlim(AGE_MAX, AGE_MIN)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("Age (Ma)")

    def boxes(items, y0, height, fontsize):
        for name, abbr, top, base in items:
            x0, x1 = max(top, AGE_MIN), min(base, AGE_MAX)
            if x1 <= x0:
                continue
            ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, height,
                                       facecolor="white", edgecolor="black",
                                       lw=0.6))
            xc = (x0 + x1) / 2
            for label in (name, abbr):
                t = ax.text(xc, y0 + height / 2, label, ha="center",
                            va="center", fontsize=fontsize)
                if _fits(ax, renderer, t, xc, x0, x1):
                    break
                t.remove()

    boxes(TIMESCALE_ERAS, 0.5, 0.5, 8.5)
    boxes(TIMESCALE_PERIODS, 0.0, 0.5, 7.5)


def make_figure(out_path, include_scotland=False):
    """Build and save the single-panel compilation figure."""
    fig, (ax, axts) = plt.subplots(
        2, 1, figsize=(12, 6.2), sharex=True,
        gridspec_kw={"height_ratios": [1, 0.12], "hspace": 0.0})
    plot_panel(ax, include_scotland=include_scotland)
    # Drop the (unlabeled) bottom ticks on the paleolatitude panel so the
    # timescale strip can sit flush beneath it without a redundant tick row.
    ax.tick_params(labelbottom=False, bottom=False)
    fig.canvas.draw()  # establish a renderer for the timescale fit test
    draw_timescale(axts, fig.canvas.get_renderer())

    fig.suptitle("Position of Duluth (lat = 46.79°N, lon = 92.10°W) implied by "
                 "Laurentia poles", fontsize=13, y=0.92)
    ax.legend(loc="lower left", fontsize=9, framealpha=0.9)

    fig.savefig(out_path + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(out_path + ".pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    os.makedirs(os.path.dirname(OUT_BASE), exist_ok=True)
    # Produced twice: the default (age-vetted) version, and a "_withScotland"
    # version that also shows the loose-age Scotland poles as open triangles.
    # The age vetting is strict (half-range < 50 Myr), so both Stoer (+/-70) and
    # Torridon (half-range exactly 50) fail it and appear only in the
    # "_withScotland" variant, drawn hollow by marker_style.
    for include_scotland in (False, True):
        suffix = "_withScotland" if include_scotland else ""
        out = f"{OUT_BASE}_{AGE_MAX}_{AGE_MIN}{suffix}"
        make_figure(out, include_scotland=include_scotland)
        print(f"Wrote {os.path.relpath(out, ROOT)}.png / .pdf")


if __name__ == "__main__":
    main()
