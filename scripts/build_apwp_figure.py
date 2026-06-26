"""Whole-path Laurentia APWP figures colored by age.

Reads ``data/nordic_summaries/nordic_summaries_combined.csv`` (the recreated-
from-site-level poles merged by ``combine_nordic_summaries.py``), rotates the
Greenland poles into the Laurentia reference frame (Roest & Srivastava, 1989),
and plots every pole with its A95 confidence ellipse colored by nominal age --
the whole-path counterpart to ``pole_tools.plot_apwp_context`` (which shows a
single pole in context). Two figures are written to ``_static/``:

- ``Laurentia_apwp_robinson``      : 540-1779 Ma on a Robinson projection.
- ``Laurentia_apwp_orthographic``  : 717-1779 Ma on an orthographic projection
  centered on the poles so they are all in view.

    python scripts/build_apwp_figure.py
"""

import os

import matplotlib

matplotlib.use("Agg")
import cartopy  # noqa: E402
import cartopy.crs as ccrs  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import pmagpy.ipmag as ipmag  # noqa: E402
import pmagpy.pmag as pmag  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
COMBINED_CSV = os.path.join(ROOT, "data", "nordic_summaries",
                            "nordic_summaries_combined.csv")
STATIC = os.path.join(ROOT, "_static")

COLORMAP = "viridis_r"

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


def make_figure(out_base, age_min, age_max, projection, central_lon=200):
    """Build and save one APWP figure; return the number of poles.

    ``projection`` is ``"robinson"`` or ``"orthographic"``. For the orthographic
    view the map is centered on the Fisher mean of the poles in the interval so
    they all fall on the visible hemisphere.
    """
    d = load_path(age_min, age_max)

    if projection == "orthographic":
        mean = ipmag.fisher_mean(dec=d["PLONG"].tolist(), inc=d["PLAT"].tolist())
        proj = ccrs.Orthographic(central_longitude=mean["dec"],
                                  central_latitude=mean["inc"])
        figsize = (9, 9)
    else:
        proj = ccrs.Robinson(central_longitude=central_lon)
        figsize = (12, 8)

    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection=proj)
    ax.set_global()
    ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="tan",
                   edgecolor="black", linewidth=0.3)
    ax.gridlines(color="gray", linewidth=0.4, linestyle=":")

    ipmag.plot_poles_colorbar(
        ax, d["PLONG"].tolist(), d["PLAT"].tolist(), d["A95"].tolist(),
        d["nominal age"].tolist(), age_min, age_max, colormap=COLORMAP,
        markersize=30, colorbar_label="pole age (Ma)")

    for _, r in d.iterrows():
        ax.text(r["PLONG"] + 2, r["PLAT"] + 2, str(int(r["nominal age"])),
                transform=ccrs.PlateCarree(), fontsize=5.5, color="0.25")

    ax.set_title(f"Laurentia apparent polar wander path, {age_max}-{age_min} Ma "
                 f"(recreated-from-site-level poles; Greenland rotated)",
                 fontsize=12)

    fig.savefig(out_base + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    plt.close(fig)
    return len(d)


FIGURES = [
    # (out_suffix, age_min, age_max, projection)
    ("Laurentia_apwp_robinson", 540, 1779, "robinson"),
    ("Laurentia_apwp_orthographic", 717, 1779, "orthographic"),
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
