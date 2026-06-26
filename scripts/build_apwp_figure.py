"""Whole-path Laurentia APWP figures colored by age.

Reads ``data/nordic_summaries/nordic_summaries_combined.csv`` (the recreated-
from-site-level poles merged by ``combine_nordic_summaries.py``), rotates the
Greenland and Scotland poles into the Laurentia reference frame (Greenland: Roest
& Srivastava, 1989; Scotland: Torsvik & Cocks, 2017), and plots every pole with
its A95 confidence ellipse colored by nominal age -- the whole-path counterpart
to ``pole_tools.plot_apwp_context`` (which shows a single pole in context). Each
terrane is drawn with its own marker (Laurentia circle, Greenland square, Scotland
triangle). Three figures are written to ``_static/``:

- ``Laurentia_apwp_robinson``      : 540-1779 Ma on a Robinson projection.
- ``Laurentia_apwp_orthographic``  : 717-1779 Ma on an orthographic projection
  centered on the poles so they are all in view.
- ``Laurentia_apwp_lambert``       : 717-1779 Ma on a Lambert azimuthal equal-area
  projection centered on and cropped to the pole cluster, so the A95 ellipses stay
  area-comparable with minimal wasted space.

    python scripts/build_apwp_figure.py
"""

import os

import matplotlib

matplotlib.use("Agg")
import cartopy  # noqa: E402
import cartopy.crs as ccrs  # noqa: E402
import cmcrameri.cm  # noqa: E402,F401  registers the 'cmc.*' colormaps
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import pmagpy.ipmag as ipmag  # noqa: E402
import pmagpy.pmag as pmag  # noqa: E402
from adjustText import adjust_text  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
COMBINED_CSV = os.path.join(ROOT, "data", "nordic_summaries",
                            "nordic_summaries_combined.csv")
STATIC = os.path.join(ROOT, "_static")

COLORMAP = "Spectral_r"

# Marker shape and size per terrane group (color still encodes age via the
# colorbar). A separate shape legend distinguishes the cratonic blocks. The list
# is ordered most-specific-first because terrane_group() matches by substring;
# the legend is displayed in LEGEND_ORDER instead. The Greenland square is given
# a smaller size because filled squares read perceptually larger than circles and
# triangles of equal point size.
TERRANE_MARKERS = [
    ("Scotland", "Scotland", "^", 30),    # rotated into Laurentia frame
    ("Greenland", "Greenland", "s", 22),  # rotated into Laurentia frame
    ("Laurentia", "Laurentia", "o", 30),  # catch-all; must be tested last
]
LEGEND_ORDER = ["Laurentia", "Greenland", "Scotland"]


def terrane_group(terrane):
    """Map a Terrane string to (label, marker, size), tested most-specific first."""
    t = str(terrane)
    for key, label, marker, size in TERRANE_MARKERS:
        if key in t:
            return label, marker, size
    return "Laurentia", "o", 30

# Euler poles [pole_lat, pole_lon, angle] used to rotate separated-terrane poles
# into the Laurentia reference frame, mirroring pole_tools.TERRANE_EULER_POLES
# (Greenland: Roest & Srivastava, 1989; Scotland: Torsvik & Cocks, 2017).
# Svalbard is excluded from the figure.
TERRANE_EULER_POLES = {
    "Laurentia-Greenland": [67.5, -118.5, -13.8],
    "Laurentia-Greenland-Nain": [67.5, -118.5, -13.8],
    "Laurentia-Scotland": [78.6, 161.9, -31.0],
}


def load_path(age_min, age_max):
    """Load the combined summaries, rotate Greenland poles, restrict to the
    age interval, sorted by age."""
    d = pd.read_csv(COMBINED_CSV)
    for c in ["PLAT", "PLONG", "A95", "nominal age"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["PLAT", "PLONG", "A95", "nominal age"]).copy()
    d = d[~d["Terrane"].astype(str).str.contains("Svalbard")]
    d = d[(d["nominal age"] >= age_min) & (d["nominal age"] <= age_max)].copy()

    # rotate separated-terrane (Greenland, Scotland) poles into the Laurentia frame
    for idx in d.index:
        euler = TERRANE_EULER_POLES.get(str(d.at[idx, "Terrane"]))
        if euler is None:
            continue
        rlat, rlon = pmag.pt_rot(euler, [d.at[idx, "PLAT"]], [d.at[idx, "PLONG"]])
        d.at[idx, "PLAT"], d.at[idx, "PLONG"] = rlat[0], rlon[0]

    return d.sort_values("nominal age").reset_index(drop=True)


def cluster_extent(proj, d, pad_deg=16):
    """Projected-coordinate bounds enclosing the poles, padded by ``pad_deg``
    (degrees, converted to metres) to leave room for the A95 ellipses and labels."""
    pts = proj.transform_points(ccrs.PlateCarree(),
                                d["PLONG"].to_numpy(), d["PLAT"].to_numpy())
    x, y = pts[:, 0], pts[:, 1]
    pad = pad_deg * 111320.0  # deg -> m, approximate
    return [x.min() - pad, x.max() + pad, y.min() - pad, y.max() + pad]


def make_figure(out_base, age_min, age_max, projection, central_lon=200):
    """Build and save one APWP figure; return the number of poles.

    ``projection`` is ``"robinson"``, ``"orthographic"``, or ``"lambert"``. The
    orthographic and lambert views are centered on the Fisher mean of the poles in
    the interval; the lambert (equal-area) view is additionally cropped to the
    pole cluster so the A95 ellipses stay area-comparable with no wasted space.
    """
    d = load_path(age_min, age_max)

    if projection in ("orthographic", "lambert"):
        mean = ipmag.fisher_mean(dec=d["PLONG"].tolist(), inc=d["PLAT"].tolist())
        if projection == "lambert":
            proj = ccrs.LambertAzimuthalEqualArea(
                central_longitude=mean["dec"], central_latitude=mean["inc"])
        else:
            proj = ccrs.Orthographic(central_longitude=mean["dec"],
                                     central_latitude=mean["inc"])
        figsize = (11, 9) if projection == "lambert" else (9, 9)
    else:
        proj = ccrs.Robinson(central_longitude=central_lon)
        figsize = (12, 8)

    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection=proj)
    if projection == "lambert":
        ax.set_extent(cluster_extent(proj, d), crs=proj)
    else:
        ax.set_global()
    ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="tan",
                   edgecolor="black", linewidth=0.3)
    ax.gridlines(color="gray", linewidth=0.4, linestyle=":")

    d["_label"] = [terrane_group(t)[0] for t in d["Terrane"]]

    # plot each terrane group with its own marker/size; draw the colorbar once
    first = True
    for _, label, marker, size in TERRANE_MARKERS:
        g = d[d["_label"] == label]
        if g.empty:
            continue
        ipmag.plot_poles_colorbar(
            ax, g["PLONG"].tolist(), g["PLAT"].tolist(), g["A95"].tolist(),
            g["nominal age"].tolist(), age_min, age_max, colormap=COLORMAP,
            marker=marker, markersize=size, colorbar=first,
            colorbar_label="pole age (Ma)")
        first = False

    # age labels, de-overlapped with adjustText. Labels are placed in the axes'
    # projected (transData) coordinate system -- not PlateCarree -- so adjustText
    # can move them in a flat plane; thin leader lines connect each displaced
    # label back to its pole.
    xs, ys, texts = [], [], []
    for _, r in d.iterrows():
        x, y = proj.transform_point(r["PLONG"], r["PLAT"], ccrs.PlateCarree())
        xs.append(x)
        ys.append(y)
        texts.append(ax.text(x, y, str(int(r["nominal age"])), fontsize=8,
                             color="0.15", ha="center", va="center",
                             bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                                       edgecolor="none", alpha=0.7)))
    adjust_text(texts, x=xs, y=ys, ax=ax,
                force_text=(0.3, 0.5), expand=(1.3, 1.5),
                arrowprops=dict(arrowstyle="-", color="0.5", lw=0.4))

    # shape legend (neutral color; pole color encodes age). Legend marker sizes
    # are scaled from the plotted sizes so the square stays proportionally smaller.
    style_for = {lab: (mk, sz) for _, lab, mk, sz in TERRANE_MARKERS}
    handles = [plt.Line2D([], [], linestyle="none", marker=style_for[lab][0],
                          color="0.4", markeredgecolor="k",
                          markersize=style_for[lab][1] * 8 / 30, label=lab)
               for lab in LEGEND_ORDER if (d["_label"] == lab).any()]
    ax.legend(handles=handles, loc="lower right", frameon=True, fontsize=12)

    fig.savefig(out_base + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    plt.close(fig)
    return len(d)


FIGURES = [
    # (out_suffix, age_min, age_max, projection)
    ("Laurentia_apwp_robinson", 540, 1779, "robinson"),
    ("Laurentia_apwp_orthographic", 717, 1779, "orthographic"),
    ("Laurentia_apwp_lambert", 717, 1779, "lambert"),
]


def main():
    os.makedirs(STATIC, exist_ok=True)
    for suffix, age_min, age_max, projection in FIGURES:
        out = os.path.join(STATIC, suffix)
        n = make_figure(out, age_min, age_max, projection)
        print(f"Wrote {os.path.relpath(out, ROOT)}.png / .pdf "
              f"({age_max}-{age_min} Ma, {projection}: {n} poles)")


if __name__ == "__main__":
    main()
