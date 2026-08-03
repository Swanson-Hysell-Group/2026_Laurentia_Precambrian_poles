"""Whole-path Laurentia APWP figures colored by age.

Reads ``data/nordic_summaries/nordic_summaries_combined.csv`` (the recreated-
from-site-level poles merged by ``combine_nordic_summaries.py``), rotates the
Greenland and Scotland poles into the Laurentia reference frame (Greenland: Roest
& Srivastava, 1989; Scotland: Torsvik & Cocks, 2017), and plots every pole with
its A95 confidence ellipse colored by nominal age -- the whole-path counterpart
to ``pole_tools.plot_apwp_context`` (which shows a single pole in context). Each
terrane is drawn with its own marker (Laurentia circle, Greenland square, Scotland
triangle). The sedimentary units are plotted at their inclination-shallowing-
corrected (Kent mean) positions from
``data/nordic_summaries/kent_poles_combined.csv``, drawn with their true Kent
95% ellipse rather than a circular A95.

Two SphereUDE fit inputs are exported, differing only in how the sedimentary
poles enter (see ``TRACKS``): ``apwp_fit_input_uncorrected.csv`` takes them as
measured, and ``apwp_fit_input_corrected.csv`` takes the inclination-
shallowing-corrected Kent means with the ellipse simplified to its equal-area
circular A95, which is the scalar SphereUDE's Fisher weighting can consume. The
figures overlay the ``corrected`` path, matching the poles they plot. Three
figures are written to ``_static/``:

- ``Laurentia_apwp_robinson``      : 540-1779 Ma on a Robinson projection.
- ``Laurentia_apwp_orthographic``  : 717-1779 Ma on an orthographic projection
  centered on the poles so they are all in view.
- ``Laurentia_apwp_lambert``       : 717-1779 Ma on a Lambert azimuthal equal-area
  projection centered on and cropped to the pole cluster, so the A95 ellipses stay
  area-comparable with minimal wasted space.

    python scripts/build_apwp_figure.py
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import cartopy  # noqa: E402
import cartopy.crs as ccrs  # noqa: E402
import cmcrameri.cm  # noqa: E402,F401  registers the 'cmc.*' colormaps
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pmagpy.ipmag as ipmag  # noqa: E402
import pmagpy.pmag as pmag  # noqa: E402
from adjustText import adjust_text  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from build_pole_map import (  # noqa: E402
    KENT_ELLIPSE_COLUMNS, apply_kent_poles, kent_dict)

ROOT = os.path.dirname(HERE)
COMBINED_CSV = os.path.join(ROOT, "data", "nordic_summaries",
                            "nordic_summaries_combined.csv")
STATIC = os.path.join(ROOT, "_static")

# --- the two APWP tracks ----------------------------------------------------
# The compilation supports two spline paths, differing only in how the
# sedimentary poles enter:
#
#   "uncorrected" = the poles as measured, i.e. detrital inclinations still
#                   shallowed by compaction.
#   "corrected"   = the inclination-shallowing-corrected (Kent mean) positions.
#
# Named by whether the correction has been applied, rather than by
# "flattened"/"unflattened", because the latter pair inverts easily in reading
# (compaction flattens; unflattening is the correction, so "unflattened" means
# the *corrected* pole). Note that ``pole_tools.make_nordic_summary``'s
# ``pole_mean_unflattened`` argument and ``unflattened_pole`` in the Lower Freda
# notebook still use that older vocabulary for the corrected pole.
#
# SphereUDE weights each pole by a single Fisher concentration
# ``kappa ~ (140/A95)**2``, so it cannot consume a Kent ellipse's two
# concentration parameters. The corrected input therefore carries the
# equal-area circular simplification of the ellipse, ``A95 = sqrt(zeta95 *
# eta95)``, which ``apply_kent_poles`` has already put in the ``A95`` column;
# any directionality in the Kent confidence region is discarded by the fit.
TRACKS = {
    "uncorrected": dict(
        kent=False,
        fit_input="apwp_fit_input_uncorrected.csv",
        path="apwp_sphereude_path_uncorrected.csv",
        label="poles as measured (no inclination-shallowing correction)"),
    "corrected": dict(
        kent=True,
        fit_input="apwp_fit_input_corrected.csv",
        path="apwp_sphereude_path_corrected.csv",
        label="inclination-shallowing-corrected (Kent mean) sedimentary poles"),
}
# The track the figures plot. The poles drawn on the maps are the Kent-corrected
# ones, so the overlaid path is the one fit to them.
FIGURE_TRACK = "corrected"


def track_path(name, key):
    """Absolute path of a track's ``fit_input`` or ``path`` CSV."""
    return os.path.join(ROOT, "data", "nordic_summaries", TRACKS[name][key])

COLORMAP = "cmc.managua"

# The ca. 1382 Ma Greenland poles (Midsommersoe, Victoria Fjord, Zig-Zag Dal) are
# still plotted but excluded from the path fit: they conflict with the rest of the
# path.
FIT_EXCLUDE_GREENLAND_AGE = 1382
# Poles whose age is uncertain by more than this half-range (Myr) are excluded
# from the fit so loosely-dated poles do not pull on the path. At 50 Myr this
# drops the two Scotland poles (Torridon, Stoer) and nothing else.
FIT_MAX_AGE_HALF_UNC = 50

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


def plot_kent_ellipse(map_axis, kent, **kwargs):
    """Draw a Kent mean pole and its 95% ellipse on the correct hemisphere.

    ``ipmag.plot_pole_ellipse`` passes ``lower`` through to
    ``pmagplotlib.plot_ell``, whose default (``lower=True``) returns the
    lower-hemisphere projection of the ellipse. For a pole at a **southern**
    latitude that projection is the antipode, so the ellipse is drawn ~170 deg
    from the pole it belongs to -- and silently, since the pole marker is
    plotted separately and still lands in the right place. Choosing ``lower``
    from the sign of the pole latitude keeps the ellipse around its own pole.
    ``pole_tools.compilation_kent_pole`` applies the same correction to the
    diagnostic map it draws in the notebooks.

    Args:
        map_axis: Cartopy axis to draw on.
        kent (dict): Kent statistics, e.g. from
            :func:`build_pole_map.kent_dict`.
        **kwargs: Passed to ``ipmag.plot_pole_ellipse``.

    Returns:
        The map axis.
    """
    return ipmag.plot_pole_ellipse(map_axis, kent,
                                   lower=kent["inc"] >= 0, **kwargs)


def _no_kent(d):
    """Add the ``_kent``/ellipse columns as empty, leaving every pole as measured.

    Keeps the ``uncorrected`` track's frame shape identical to the ``corrected``
    one so the rotation and plotting code below is shared.
    """
    d = d.copy()
    d["_kent"] = False
    d["_kent_method"] = ""
    for col in KENT_ELLIPSE_COLUMNS:
        d["_" + col] = np.nan
    return d


def load_path(age_min, age_max, kent=True):
    """Load the combined summaries, rotate Greenland poles, restrict to the
    age interval, sorted by age.

    Args:
        age_min, age_max (float): Age interval to keep, in Ma.
        kent (bool): With ``True`` (the ``corrected`` track) the sedimentary
            units are substituted with their inclination-shallowing-corrected
            Kent mean poles, the ``A95`` column carrying the equal-area circular
            simplification of the Kent ellipse. With ``False`` (the
            ``uncorrected`` track) every pole is left as measured.

    Returns:
        pd.DataFrame: Poles sorted by nominal age, rotated into the Laurentia
        frame.
    """
    d = pd.read_csv(COMBINED_CSV)
    for c in ["PLAT", "PLONG", "A95", "nominal age"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = apply_kent_poles(d, verbose=kent) if kent else _no_kent(d)
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
        # the Kent ellipse axes are directions on the sphere, so they rotate
        # with the pole; the semi-angles (Zeta, Eta) are unchanged by a rigid
        # rotation
        if not d.at[idx, "_kent"]:
            continue
        for dec_col, inc_col in (("_Zdec", "_Zinc"), ("_Edec", "_Einc")):
            alat, alon = pmag.pt_rot(euler, [d.at[idx, inc_col]],
                                     [d.at[idx, dec_col]])
            d.at[idx, inc_col], d.at[idx, dec_col] = alat[0], alon[0]

    return d.sort_values("nominal age").reset_index(drop=True)


def export_fit_input(age_min, age_max, track=FIGURE_TRACK, out_csv=None):
    """Write the rotated, filtered poles used as input to the SphereUDE path fit.

    Starts from ``load_path`` (Greenland and Scotland already rotated into the
    Laurentia frame, Svalbard dropped) and additionally drops the ca. 1382 Ma
    Greenland poles, which conflict with the rest of the path, and any pole whose
    age half-range exceeds ``FIT_MAX_AGE_HALF_UNC``.

    The four columns written (``age``, ``plat``, ``plon``, ``a95``) are the same
    for both tracks; only the sedimentary rows differ. For the ``corrected``
    track ``a95`` is the equal-area circular simplification of the Kent ellipse,
    which is the scalar SphereUDE's Fisher weighting can consume.

    Args:
        age_min, age_max (float): Age interval to fit.
        track (str): Key of :data:`TRACKS` -- ``'uncorrected'`` or
            ``'corrected'``.
        out_csv (str or None): Output path; defaults to the track's own file.

    Returns:
        int: Number of poles written.
    """
    if out_csv is None:
        out_csv = track_path(track, "fit_input")
    d = load_path(age_min, age_max, kent=TRACKS[track]["kent"])
    # exclude the ca. 1382 Ma Greenland poles (conflict with the rest of the path)
    drop = (d["Terrane"].astype(str).str.contains("Greenland") &
            (d["nominal age"] == FIT_EXCLUDE_GREENLAND_AGE))
    # exclude loosely-dated poles (age half-range > FIT_MAX_AGE_HALF_UNC Myr)
    half_unc = (pd.to_numeric(d["himagage"], errors="coerce") -
                pd.to_numeric(d["lomagage"], errors="coerce")) / 2
    drop = drop | (half_unc > FIT_MAX_AGE_HALF_UNC)
    d = d.loc[~drop]
    d[["nominal age", "PLAT", "PLONG", "A95"]].rename(
        columns={"nominal age": "age", "PLAT": "plat",
                 "PLONG": "plon", "A95": "a95"}).to_csv(out_csv, index=False)
    return len(d)


def load_spline_path(track=FIGURE_TRACK):
    """Load a SphereUDE fitted path (age, lat, lon), or ``None`` if absent.

    Each track is read only from its own file. There is deliberately no fallback
    to the pre-track name ``apwp_sphereude_path.csv``: that file is not
    self-describing, so a fallback silently mislabels whichever pole set it
    happens to hold. A missing path is not an error -- the figures are simply
    drawn without the spline (see ``scripts/sphereude/README.md``).
    """
    candidate = track_path(track, "path")
    if os.path.exists(candidate):
        return pd.read_csv(candidate)
    return None


def plot_age_graded_path(ax, path, vmin, vmax, lw=3.5, zorder=3):
    """Draw the fitted path as a thick line whose color progresses along the age
    colormap. Segments are split where the longitude wraps so the line does not
    streak across the map."""
    from matplotlib.collections import LineCollection

    lon = path["lon"].to_numpy()
    lat = path["lat"].to_numpy()
    age = path["age"].to_numpy()
    pts = ax.projection.transform_points(ccrs.PlateCarree(), lon, lat)[:, :2]
    # break segments where the projected step is implausibly large (wrap-around)
    step = np.hypot(np.diff(pts[:, 0]), np.diff(pts[:, 1]))
    ok = step < 5 * np.median(step)
    # Overlapping three-point segments rather than abutting two-point ones. Each
    # segment of a LineCollection is stroked separately, so where two merely abut
    # their antialiased edges composite against the background instead of against
    # each other and leave a hairline seam at every join. That is invisible in
    # raster output at figure dpi but reads as a dashed line in the vector PDF.
    # Overlapping by a full step means each segment paints over its predecessor's
    # join, so the drawn path is continuous.
    segs = np.stack([pts[:-2], pts[1:-1], pts[2:]], axis=1)
    keep = ok[:-1] & ok[1:]          # drop segments touching a wrap-around step
    lc = LineCollection(segs[keep], cmap=COLORMAP, zorder=zorder,
                        linewidth=lw, capstyle="round", joinstyle="round")
    lc.set_array(age[1:-1][keep])
    lc.set_clim(vmin, vmax)
    ax.add_collection(lc)


def cluster_extent(proj, d, pad_deg=16):
    """Projected-coordinate bounds enclosing the poles, padded by ``pad_deg``
    (degrees, converted to metres) to leave room for the A95 ellipses and labels."""
    pts = proj.transform_points(ccrs.PlateCarree(),
                                d["PLONG"].to_numpy(), d["PLAT"].to_numpy())
    x, y = pts[:, 0], pts[:, 1]
    pad = pad_deg * 111320.0  # deg -> m, approximate
    return [x.min() - pad, x.max() + pad, y.min() - pad, y.max() + pad]


def add_named_path_annotations(ax, proj):
    """Annotate named late Mesoproterozoic segments of the Laurentia APWP.

    Places three labels -- the Logan Loop (ca. 1144-1108 Ma hairpin), the
    Keweenawan Track (ca. 1108-1080 Ma descending limb), and the Grenville Loop
    (ca. 1080-990 Ma broad loop) -- in open areas of the Lambert figure, each
    connected to its path segment by a subtle leader line. Target and label
    positions are given in longitude-latitude and transformed into the Lambert
    projection. Returns the list of annotation artists so they can be passed to
    ``adjust_text(objects=...)`` and thereby repel the numerical age labels.
    """
    annotations = [
        {
            "text": "Logan Loop",
            "target": (229, 64),
            "label": (229, 64),
            "rotation": 0,
        },
        {
            "text": "Grenville Loop",
            "target": (195, 5),
            "label": (162, -23),
            "rotation": 0,
        },
    ]

    artists = []
    for item in annotations:
        target_x, target_y = proj.transform_point(
            item["target"][0], item["target"][1], ccrs.PlateCarree())
        label_x, label_y = proj.transform_point(
            item["label"][0], item["label"][1], ccrs.PlateCarree())
        artist = ax.annotate(
            item["text"],
            xy=(target_x, target_y),
            xytext=(label_x, label_y),
            xycoords="data",
            textcoords="data",
            fontsize=10,
            fontweight="semibold",
            color="0.15",
            ha="center",
            va="center",
            rotation=item["rotation"],
            rotation_mode="anchor",
            bbox=dict(
                boxstyle="round,pad=0.2",
                facecolor="white",
                edgecolor="none",
                alpha=0.85,
            ),
            zorder=8,
        )
        artists.append(artist)
    return artists


def make_figure(out_base, age_min, age_max, projection, central_lon=200):
    """Build and save one APWP figure; return the number of poles.

    ``projection`` is ``"robinson"``, ``"orthographic"``, or ``"lambert"``. The
    orthographic and lambert views are centered on the Fisher mean of the poles in
    the interval; the lambert (equal-area) view is additionally cropped to the
    pole cluster so the A95 ellipses stay area-comparable with no wasted space.
    """
    d = load_path(age_min, age_max, kent=TRACKS[FIGURE_TRACK]["kent"])

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
                   edgecolor="none")
    if projection == "lambert":
        # label longitude along the bottom and latitude along the left edge.
        # Longitudes are shown in 0-360 degrees east (the convention poles are
        # reported in); latitude gridlines are placed every 15 degrees.
        from matplotlib.ticker import FuncFormatter, MultipleLocator
        gl = ax.gridlines(color="gray", linewidth=0.4, linestyle=":",
                          draw_labels={"bottom": "x", "left": "y"})
        gl.ylocator = MultipleLocator(15)
        gl.xformatter = FuncFormatter(lambda lon, _: f"{int(round(lon % 360))}°E")
        gl.xlabel_style = {"size": 9, "color": "0.3"}
        gl.ylabel_style = {"size": 9, "color": "0.3"}
    else:
        ax.gridlines(color="gray", linewidth=0.4, linestyle=":")

    # age-graded SphereUDE path (drawn under the markers), if it has been fitted
    path = load_spline_path(FIGURE_TRACK)
    if path is not None:
        plot_age_graded_path(ax, path, age_min, age_max)

    # named APWP segments (Lambert view only), created before the poles and age
    # labels so adjustText can be told to route the age labels around them
    named_annotations = []
    if projection == "lambert":
        named_annotations = add_named_path_annotations(ax, proj)

    d["_label"] = [terrane_group(t)[0] for t in d["Terrane"]]

    # Plot each terrane group with its own marker/size; draw the colorbar once.
    # The sedimentary units carry a Kent confidence ellipse rather than a
    # circular A95, so they are drawn with ipmag.plot_pole_ellipse (which takes
    # the ellipse's own axes) instead of the circular plot_poles_colorbar; the
    # color still encodes age, taken from the same normalization.
    age_color = plt.cm.ScalarMappable(
        cmap=COLORMAP, norm=plt.Normalize(vmin=age_min, vmax=age_max))
    first = True
    for _, label, marker, size in TERRANE_MARKERS:
        g = d[d["_label"] == label]
        if g.empty:
            continue
        circular = g[~g["_kent"]]
        if not circular.empty:
            ipmag.plot_poles_colorbar(
                ax, circular["PLONG"].tolist(), circular["PLAT"].tolist(),
                circular["A95"].tolist(), circular["nominal age"].tolist(),
                age_min, age_max, colormap=COLORMAP,
                marker=marker, markersize=size, colorbar=first,
                colorbar_label="pole age (Ma)")
            first = False
        for _, r in g[g["_kent"]].iterrows():
            plot_kent_ellipse(
                ax, kent_dict(r),
                color=age_color.to_rgba(float(r["nominal age"])),
                edgecolor="k", marker=marker, markersize=size)

    # age labels, de-overlapped with adjustText. Labels are placed in the axes'
    # projected (transData) coordinate system -- not PlateCarree -- so adjustText
    # can move them in a flat plane; thin leader lines connect each displaced
    # label back to its pole.
    xs, ys = [], []
    for _, r in d.iterrows():
        x, y = proj.transform_point(r["PLONG"], r["PLAT"], ccrs.PlateCarree())
        xs.append(x)
        ys.append(y)
    # nudge each label slightly to the left of its pole so the default resting
    # place clears the marker; adjustText still anchors leader lines to xs/ys.
    dx = 0.015 * (max(xs) - min(xs))
    texts = [ax.text(x - dx, y, str(int(r["nominal age"])), fontsize=8,
                     color="0.15", ha="center", va="center", zorder=102,
                     bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                               edgecolor="none", alpha=0.85))
             for (x, y), (_, r) in zip(zip(xs, ys), d.iterrows())]
    adjust_text(texts, x=xs, y=ys, ax=ax,
                objects=named_annotations if named_annotations else None,
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

    # Export the poles for the SphereUDE path fit (717-1779 Ma, the cluster
    # range), once per track. The two inputs differ only in the sedimentary
    # rows; everything else -- rotation, exclusions, age vetting -- is shared.
    for name, spec in TRACKS.items():
        n_fit = export_fit_input(717, 1779, track=name)
        print(f"Wrote {os.path.relpath(track_path(name, 'fit_input'), ROOT)} "
              f"({n_fit} poles, {name}: {spec['label']})")
        if load_spline_path(name) is None:
            print(f"  (no fitted path for this track yet; fit it to "
                  f"{spec['fit_input']} and write "
                  f"{spec['path']})")
    print(f"Figures overlay the '{FIGURE_TRACK}' path, matching the poles they "
          f"plot")

    for suffix, age_min, age_max, projection in FIGURES:
        out = os.path.join(STATIC, suffix)
        n = make_figure(out, age_min, age_max, projection)
        print(f"Wrote {os.path.relpath(out, ROOT)}.png / .pdf "
              f"({age_max}-{age_min} Ma, {projection}: {n} poles)")


if __name__ == "__main__":
    main()
