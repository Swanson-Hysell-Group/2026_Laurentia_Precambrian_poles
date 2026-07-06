"""Paleomagnetic reconstructions of Laurentia through time on a paleolatitude ladder.

Each reconstruction is a *pure paleomagnetic* reconstruction driven by the
SphereUDE spline APWP (``data/nordic_summaries/apwp_sphereude_path.csv``, fit by
``scripts/sphereude/fit_apwp_spline.jl``): for a given reconstruction age the
spline pole (plat, plon) is sampled and Laurentia is rotated so that pole sits at
the spin axis. Paleolatitude and azimuthal orientation are therefore constrained
by the path; paleolongitude is not.

To keep distortion low even when the continent sits at high paleolatitude (e.g.
the ca. 1150 Ma slice), each reconstruction is drawn in its own Lambert azimuthal
equal-area projection centered on that continent (an equal-area, minimally
distorted "globe" view, like a Mollweide centered on the continent). The panels
are then stacked on a shared vertical paleolatitude axis -- each continent's
centroid sits at its true paleolatitude -- and spaced left-to-right from oldest to
youngest. Horizontal position is arbitrary (paleolongitude is unconstrained); the
sequence simply reads as a march through time.

The reconstructed geometry is the Whitmeyer & Karlstrom (2007) basement-province
outline, drawn with the same age-category colors used on the compilation map.
Because Laurentia grows by accretion, each unit is only drawn at reconstruction
ages at or younger than its ``appearance_age`` (from
``build_province_aged_geojson.py``: the first-formation age for accreted crust,
or the inherited host age for intraplate/reworking units such as the Midcontinent
Rift), so older slices show a smaller proto-Laurentia that builds outward.

Output (written to ``_static/``):

- ``Laurentia_reconstructions_ladder.pdf`` / ``.png``

Run in an env with pmagpy/shapely (e.g. ``oman_dikes``)::

    python scripts/build_reconstruction_figure.py
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pmagpy.pmag as pmag  # noqa: E402
import shapely.geometry as sgeom  # noqa: E402
from matplotlib.patches import Patch, PathPatch  # noqa: E402
from matplotlib.path import Path  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

from build_province_geojson import CATEGORY_COLOR  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SPHEREUDE_PATH_CSV = os.path.join(ROOT, "data", "nordic_summaries",
                                  "apwp_sphereude_path.csv")
# Age-resolved provinces (one feature per unit with an ``appearance_age``), built
# by build_province_aged_geojson.py so the continent grows by accretion.
PROVINCES_GEOJSON = os.path.join(ROOT, "data", "geologic_provinces",
                                 "Whitmeyer2007_provinces_aged.geojson")
STATIC = os.path.join(ROOT, "_static")

# Reconstruction ages (Ma), oldest first so the ladder reads oldest -> youngest
# left to right. All lie within the 719-1779 Ma span of the spline path.
RECON_AGES = list(range(1750, 749, -100))

# Horizontal pitch between successive panels, as a fraction of the widest
# projected continent. Because the panels drift in paleolatitude they can be
# packed closer than their full width without colliding; lower this to tighten
# the figure, raise it toward 1 to spread the panels out.
PITCH_FRACTION = 0.66
# Extra room (deg) left of the first panel and right of the last, beyond the
# panels' own extent; the larger left pad keeps the oldest panel off the axis.
LEFT_PAD = 12.0
RIGHT_PAD = 2.0
# White backing silhouette drawn behind each panel so overlapping continents stay
# legible. HALO_BUFFER is the tight outer margin (projection degrees) the white
# extends beyond the outline. CLOSE_BUFFER is a morphological closing applied
# first (buffer out then back in) so detached pieces within ~2x its value -- e.g.
# Greenland across the Nares Strait gap -- merge into one silhouette with the
# mainland instead of getting their own enclosure.
HALO_BUFFER = 1.0
CLOSE_BUFFER = 4.5
# Drop projected polygons smaller than this (deg^2). The merged province geometry
# carries near-degenerate slivers whose buffered halo inflates into stray
# circles/dangles on the margin; real provinces are orders of magnitude larger.
MIN_POLY_AREA = 1.0
# Paleolatitude graticule lines (deg); the equator is emphasized.
LAT_GRID = [-60, -30, 0, 30, 60]


def reconstruction_euler(plat, plon):
    """Euler pole [lat, lon, angle] that sends paleopole (plat, plon) to the spin
    axis -- the pure paleomagnetic reconstruction rotation. The resulting
    paleolongitude is arbitrary (the panel is re-centered on the continent)."""
    return [0.0, plon - 90.0, 90.0 - plat]


def sample_pole(path, age):
    """Interpolate the spline (plat, plon) at a reconstruction age. Longitude is
    interpolated through Cartesian components to avoid wrap artifacts."""
    plat = np.interp(age, path["age"], path["lat"])
    lon_rad = np.radians(path["lon"].to_numpy())
    cos_i = np.interp(age, path["age"], np.cos(lon_rad))
    sin_i = np.interp(age, path["age"], np.sin(lon_rad))
    plon = np.degrees(np.arctan2(sin_i, cos_i)) % 360.0
    return float(plat), float(plon)


def spherical_centroid(lats, lons):
    """Centroid (lat0, lon0) of points as the direction of their mean unit vector."""
    la, lo = np.radians(lats), np.radians(lons)
    x = np.mean(np.cos(la) * np.cos(lo))
    y = np.mean(np.cos(la) * np.sin(lo))
    z = np.mean(np.sin(la))
    return (np.degrees(np.arctan2(z, np.hypot(x, y))),
            np.degrees(np.arctan2(y, x)))


def laea(lats, lons, lat0, lon0):
    """Lambert azimuthal equal-area projection (deg-equivalent units) centered on
    (lat0, lon0). At the center the scale is ~1 deg of arc per output unit, so the
    panel can be stacked directly against a degree-valued paleolatitude axis."""
    la, lo = np.radians(lats), np.radians(lons)
    la0, lo0 = np.radians(lat0), np.radians(lon0)
    dlon = lo - lo0
    denom = 1.0 + np.sin(la0) * np.sin(la) + np.cos(la0) * np.cos(la) * np.cos(dlon)
    k = np.sqrt(2.0 / np.clip(denom, 1e-12, None))
    x = k * np.cos(la) * np.sin(dlon)
    y = k * (np.cos(la0) * np.sin(la) - np.sin(la0) * np.cos(la) * np.cos(dlon))
    scale = np.degrees(1.0)  # 180/pi: output units ~ degrees of arc near center
    return scale * x, scale * y


def iter_rings(geom):
    """Yield (exterior, [holes]) coordinate arrays for each polygon in a (Multi)Polygon."""
    polys = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
    for poly in polys:
        yield (np.asarray(poly.exterior.coords),
               [np.asarray(r.coords) for r in poly.interiors])


def rotate_lonlat(coords, euler):
    """Rotate an (N,2) array of (lon, lat) by the Euler pole; return (rlon, rlat)."""
    rlat, rlon = pmag.pt_rot(euler, list(coords[:, 1]), list(coords[:, 0]))
    return np.asarray(rlon), np.asarray(rlat)


def main():
    path = pd.read_csv(SPHEREUDE_PATH_CSV)
    with open(PROVINCES_GEOJSON) as fh:
        provinces = json.load(fh)["features"]
    # Category draw order (Archean under, younger over) and colors come from the
    # shared definition so they match the static compilation map.
    cat_order = {c: i for i, c in enumerate(CATEGORY_COLOR)}
    cat_color = dict(CATEGORY_COLOR)

    # --- Pass 1: reconstruct each slice, store rings as (rlon, rlat), and find the
    # continent centroid (the panel's projection center and vertical position). ---
    slices = []
    for age in RECON_AGES:
        plat, plon = sample_pole(path, age)
        euler = reconstruction_euler(plat, plon)
        rings, all_lon, all_lat = [], [], []
        for feat in provinces:
            props = feat["properties"]
            if props["appearance_age"] < age:
                continue    # host crust not yet accreted at this reconstruction age
            cat = props["category"]
            # Reworking/intraplate units read as their host crust until emplaced,
            # then switch to their own color (e.g. Midcontinent Rift after ~1.2 Ga).
            if props.get("reworking") and age > props["magmatic_age"]:
                cat = props["host_category"]
            geom = sgeom.shape(feat["geometry"])
            for ext, holes in iter_rings(geom):
                ext_ll = rotate_lonlat(ext, euler)
                hole_ll = [rotate_lonlat(h, euler) for h in holes]
                rings.append((cat, ext_ll, hole_ll))
                all_lon.append(ext_ll[0])
                all_lat.append(ext_ll[1])
        all_lon = np.concatenate(all_lon)
        all_lat = np.concatenate(all_lat)
        lat0, lon0 = spherical_centroid(all_lat, all_lon)
        px, py = laea(all_lat, all_lon, lat0, lon0)
        slices.append({
            "age": age, "plat": plat, "plon": plon,
            "rings": rings, "lat0": lat0, "lon0": lon0,
            "width": px.max() - px.min(),
            "x_min": float(px.min()), "x_max": float(px.max()),
            "y_lo": lat0 + py.min(), "y_hi": lat0 + py.max(),
        })

    pitch = PITCH_FRACTION * max(s["width"] for s in slices)
    x_centers = [i * pitch for i in range(len(slices))]

    # --- Pass 2: plot each panel as an equal-area patch on the paleolatitude axis ---
    x_lo = x_centers[0] + slices[0]["x_min"] - LEFT_PAD
    x_hi = x_centers[-1] + slices[-1]["x_max"] + RIGHT_PAD
    y_lo = min(s["y_lo"] for s in slices) - 6
    y_hi = max(s["y_hi"] for s in slices) + 12
    # Size the figure to the data aspect so equal-area panels are not squashed.
    aspect = (y_hi - y_lo) / (x_hi - x_lo)
    fig_w = 1.55 * len(slices) + 1.5
    fig, ax = plt.subplots(figsize=(fig_w, fig_w * aspect + 1.6))

    for lat in LAT_GRID:
        if not (y_lo < lat < y_hi):
            continue
        ax.axhline(lat, color="0.85", lw=(1.4 if lat == 0 else 0.8), zorder=0)
        # Centered on the line with a white box behind, drawn on top of everything.
        ax.text(x_lo + 1.5, lat, f"{lat}°", va="center", ha="left",
                fontsize=15, color="black", zorder=1000,
                bbox=dict(facecolor="white", edgecolor="none",
                          boxstyle="round,pad=0.15"))

    # zorder grows per panel so each white backing sits above the previous panel
    # (younger slices, drawn later and to the right, end up on top).
    for k, (s, xc) in enumerate(zip(slices, x_centers)):
        lat0, lon0 = s["lat0"], s["lon0"]
        base_z = 10 * k

        # Project rings to panel coordinates, dropping near-degenerate slivers
        # (by exterior area) so they neither inflate the halo nor add specks.
        kept = []
        ext_polys = []
        for cat, ext_ll, holes in s["rings"]:
            px, py = laea(ext_ll[1], ext_ll[0], lat0, lon0)
            poly = sgeom.Polygon(np.column_stack([xc + px, lat0 + py]))
            if not poly.is_valid:
                poly = poly.buffer(0)
            if poly.area < MIN_POLY_AREA:
                continue
            ext_polys.append(poly)
            kept.append((cat, ext_ll, holes))

        # White silhouette: union of the kept outlines, closed to bridge the
        # Greenland gap, then buffered tight, drawn behind the colored provinces.
        halo = (unary_union(ext_polys)
                .buffer(CLOSE_BUFFER).buffer(-CLOSE_BUFFER)  # close Greenland gap
                .buffer(HALO_BUFFER))
        for geom in (halo.geoms if halo.geom_type == "MultiPolygon" else [halo]):
            hx, hy = geom.exterior.xy
            # Prominent outline around the whole block; the white fill sits behind
            # the colored provinces, but draw the outline on top of them so it
            # cleanly frames each panel.
            ax.fill(hx, hy, facecolor="white", edgecolor="none", zorder=base_z)
            ax.plot(hx, hy, color="0.2", lw=1.2, zorder=base_z + 8)

        for cat, ext_ll, holes in kept:
            verts, codes = [], []
            for lon, lat in [ext_ll] + holes:
                px, py = laea(lat, lon, lat0, lon0)
                xy = np.column_stack([xc + px, lat0 + py])
                verts.append(xy)
                ring_codes = np.full(len(xy), Path.LINETO)
                ring_codes[0], ring_codes[-1] = Path.MOVETO, Path.CLOSEPOLY
                codes.append(ring_codes)
            patch = PathPatch(
                Path(np.concatenate(verts), np.concatenate(codes)),
                facecolor=cat_color[cat], edgecolor="0.65", lw=0.15,
                zorder=base_z + 1 + cat_order[cat])
            ax.add_patch(patch)
        ax.text(xc, s["y_hi"] + 4, f"{s['age']} Ma", ha="center", va="bottom",
                fontsize=15, fontweight="bold", zorder=base_z + 9)

    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks(LAT_GRID)
    ax.set_yticklabels([])
    ax.tick_params(length=0)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_ylabel("paleolatitude", fontsize=19)

    # The more minor, areally small categories are still drawn but omitted from
    # the legend to keep it compact.
    legend_skip = {"Belt-Purcell Supergroup", "Neoprot–Cambrian rift"}
    shown = sorted({cat for s in slices for cat, *_ in s["rings"]} - legend_skip,
                   key=lambda c: cat_order[c])
    handles = [Patch(facecolor=cat_color[c], edgecolor="0.3", label=c)
               for c in shown]
    # Legend centered below the panel field. The projection, paleolongitude
    # caveat, and data source belong in the figure caption.
    ax.legend(handles=handles, loc="upper center", ncol=len(handles), fontsize=15,
              frameon=False, bbox_to_anchor=(0.5, -0.01))

    os.makedirs(STATIC, exist_ok=True)
    for ext in ("pdf", "png"):
        out = os.path.join(STATIC, f"Laurentia_reconstructions_ladder.{ext}")
        fig.savefig(out, bbox_inches="tight", dpi=300)
        print(f"Wrote {os.path.relpath(out, ROOT)}")


if __name__ == "__main__":
    main()
