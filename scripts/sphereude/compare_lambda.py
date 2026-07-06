"""Compare SphereUDE APWP fits across the smoothness hyperparameter lambda1.

Reads a set of fitted paths produced by ``fit_apwp_spline.jl`` at different
``SU_LAMBDA1`` values (stored in ``lambda_sweep/elbow_<lam>.csv``) and produces
two diagnostic figures written to ``_static/``:

- ``Laurentia_apwp_lambda_comparison`` : one map panel per lambda1, the age-graded
  path overlaid on the poles (Lambert equal-area, as in the main figure).
- ``Laurentia_apwp_lambda_lcurve``     : the L-curve -- path roughness against
  empirical misfit -- used to locate the elbow (the balance between fitting the
  poles and keeping the path smooth; Hansen, 2001).

For each lambda1 the printed table also reports the peak angular rate (deg/Myr)
and total arc length of the path.

To (re)generate the input paths (requires Julia; see README.md), run the fit at
each lambda1 into ``lambda_sweep/``, e.g.::

    for L in 3e4 5e4 7e4 1e5 2e5 3e5; do
      SU_LAMBDA1=$L SU_OMEGADEG=2.5 SU_NITER=1500 \
        SU_OUT=scripts/sphereude/lambda_sweep/elbow_$L.csv \
        julia +lts --project=scripts/sphereude scripts/sphereude/fit_apwp_spline.jl
    done

Then, from the repository root::

    python scripts/sphereude/compare_lambda.py
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import cartopy  # noqa: E402
import cartopy.crs as ccrs  # noqa: E402
import pmagpy.ipmag as ipmag  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SWEEP_DIR = os.path.join(HERE, "lambda_sweep")
STATIC = os.path.join(ROOT, "_static")
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import build_apwp_figure as b  # noqa: E402  reuse styling + helpers

AGE_MIN, AGE_MAX = 717, 1779
# lambda1 values present in lambda_sweep/, ordered low (faithful) -> high (smooth)
LAMBDAS = [("3e4", 3e4), ("5e4", 5e4), ("7e4", 7e4),
           ("1e5", 1e5), ("2e5", 2e5), ("3e5", 3e5)]


def to_xyz(lat, lon):
    """Latitude/longitude (degrees) to Cartesian unit vectors (..., 3)."""
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo),
                     np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def path_stats(path, pole_xyz, kappa, pole_age):
    """Return (misfit, roughness, peak_rate, arc_length) for one fitted path.

    ``misfit`` is the kappa-weighted von Mises-Fisher data misfit
    Sum kappa*(1 - cos theta) between each pole and the path point at its age;
    ``roughness`` is Sum ||d^2 x||^2 over the path (a bending-energy proxy);
    ``peak_rate`` is the maximum great-circle step per Myr (deg/Myr); and
    ``arc_length`` is the total path length (deg).
    """
    px = to_xyz(path["lat"].to_numpy(), path["lon"].to_numpy())
    idx = np.clip(np.searchsorted(path["age"].to_numpy(), pole_age), 0, len(px) - 1)
    cos = np.clip((pole_xyz * px[idx]).sum(axis=1), -1, 1)
    misfit = (kappa * (1 - cos)).sum()
    d2 = px[2:] - 2 * px[1:-1] + px[:-2]
    rough = (d2 ** 2).sum()
    dot = np.clip((px[:-1] * px[1:]).sum(axis=1), -1, 1)
    step = np.degrees(np.arccos(dot))
    peak = (step / np.diff(path["age"].to_numpy())).max()
    return misfit, rough, peak, step.sum()


def main():
    os.makedirs(STATIC, exist_ok=True)
    fit = pd.read_csv(b.FIT_INPUT_CSV)
    fit["kappa"] = (140.0 / fit["a95"]) ** 2
    pole_xyz = to_xyz(fit["plat"].to_numpy(), fit["plon"].to_numpy())

    paths, stats = {}, []
    for tag, lam in LAMBDAS:
        p = pd.read_csv(os.path.join(SWEEP_DIR, f"elbow_{tag}.csv"))
        paths[tag] = p
        m, r, peak, length = path_stats(p, pole_xyz, fit["kappa"].to_numpy(),
                                        fit["age"].to_numpy())
        stats.append(dict(tag=tag, lam=lam, misfit=m, rough=r,
                          peak=peak, length=length))
    t = pd.DataFrame(stats)
    print(t.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    # --- L-curve ---
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.plot(t["misfit"], t["rough"], "o-", color="0.3")
    for _, r in t.iterrows():
        ax.annotate(f"λ1={r['lam']:.0e}", (r["misfit"], r["rough"]),
                    textcoords="offset points", xytext=(6, 4), fontsize=10)
    ax.set_xlabel(r"empirical misfit  $\Sigma\,\kappa\,(1-\cos\theta)$")
    ax.set_ylabel(r"path roughness  $\Sigma\,\|d^2x\|^2$")
    ax.set_title("APWP spline L-curve over λ1 (ωmax = 2.5°/Myr)")
    fig.tight_layout()
    fig.savefig(os.path.join(STATIC, "Laurentia_apwp_lambda_lcurve.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- path comparison grid ---
    d = b.load_path(AGE_MIN, AGE_MAX)
    d["_label"] = [b.terrane_group(x)[0] for x in d["Terrane"]]
    mean = ipmag.fisher_mean(dec=d["PLONG"].tolist(), inc=d["PLAT"].tolist())
    proj = ccrs.LambertAzimuthalEqualArea(mean["dec"], mean["inc"])

    fig = plt.figure(figsize=(21, 13))
    for i, (tag, lam) in enumerate(LAMBDAS, 1):
        ax = fig.add_subplot(2, 3, i, projection=proj)
        ax.set_extent(b.cluster_extent(proj, d), crs=proj)
        ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="tan",
                       edgecolor="black", linewidth=0.3)
        ax.gridlines(color="gray", linewidth=0.4, linestyle=":")
        b.plot_age_graded_path(ax, paths[tag], AGE_MIN, AGE_MAX, lw=3.0)
        first = True
        for _, lab, mk, sz in b.TERRANE_MARKERS:
            g = d[d["_label"] == lab]
            if g.empty:
                continue
            ipmag.plot_poles_colorbar(
                ax, g["PLONG"].tolist(), g["PLAT"].tolist(), g["A95"].tolist(),
                g["nominal age"].tolist(), AGE_MIN, AGE_MAX, colormap=b.COLORMAP,
                marker=mk, markersize=sz, colorbar=first,
                colorbar_label="pole age (Ma)")
            first = False
        peak = t.loc[t["tag"] == tag, "peak"].iloc[0]
        ax.set_title(f"λ1 = {lam:.0e}    peak {peak:.2f}°/Myr", fontsize=14)
    fig.tight_layout()
    fig.savefig(os.path.join(STATIC, "Laurentia_apwp_lambda_comparison.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Wrote _static/Laurentia_apwp_lambda_lcurve.png and "
          "_static/Laurentia_apwp_lambda_comparison.png")


if __name__ == "__main__":
    main()
