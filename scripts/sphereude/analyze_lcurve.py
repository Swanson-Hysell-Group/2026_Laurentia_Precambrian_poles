"""Assemble the L-curve from the lambda1 sweep: data misfit vs path roughness,
plus peak rate, for each lambda1. Find the elbow."""
import sys, os
sys.path.insert(0, "scripts")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.environ["SCRATCH"]
LAMBDAS = ["1e4", "3e4", "1e5", "3e5", "1e6", "3e6", "1e7"]
D = "data/nordic_summaries"

fitdf = pd.read_csv(f"{D}/apwp_fit_input.csv")
fitdf["kappa"] = (140.0 / fitdf["a95"]) ** 2


def to_xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo),
                     np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


pole_xyz = to_xyz(fitdf["plat"].values, fitdf["plon"].values)

rows = []
for L in LAMBDAS:
    f = f"{D}/lcurve_{L}.csv"
    if not os.path.exists(f):
        print(L, "MISSING")
        continue
    p = pd.read_csv(f)
    path_xyz = to_xyz(p["lat"].values, p["lon"].values)

    # kappa-weighted empirical misfit: sum kappa*(1-cos theta) at each pole's age
    idx = np.searchsorted(p["age"].values, fitdf["age"].values)
    idx = np.clip(idx, 0, len(p) - 1)
    cos = np.clip((pole_xyz * path_xyz[idx]).sum(axis=1), -1, 1)
    misfit = (fitdf["kappa"].values * (1 - cos)).sum()
    mean_dev = np.degrees(np.arccos(cos)).mean()

    # roughness proxy: sum of squared second differences of the 3D path
    d2 = path_xyz[2:] - 2 * path_xyz[1:-1] + path_xyz[:-2]
    rough = (d2 ** 2).sum()

    # peak rate (deg/Myr) and total arc length (deg)
    dot = np.clip((path_xyz[:-1] * path_xyz[1:]).sum(axis=1), -1, 1)
    step = np.degrees(np.arccos(dot))
    dt = np.diff(p["age"].values)
    peak = (step / dt).max()
    length = step.sum()

    rows.append(dict(lam=float(L), misfit=misfit, mean_dev=mean_dev,
                     rough=rough, peak=peak, length=length))

t = pd.DataFrame(rows)
print(t.to_string(index=False,
      float_format=lambda x: f"{x:.4g}"))

# L-curve plot (log-log): roughness vs misfit, labeled by lambda1
fig, ax = plt.subplots(figsize=(7, 6))
ax.loglog(t["misfit"], t["rough"], "o-", color="0.3")
for _, r in t.iterrows():
    ax.annotate(f"{r['lam']:.0e}", (r["misfit"], r["rough"]),
                textcoords="offset points", xytext=(6, 4), fontsize=9)
ax.set_xlabel("empirical misfit  Σκ(1−cosθ)")
ax.set_ylabel("path roughness  Σ‖d²x‖²")
ax.set_title("L-curve over λ1 (ωmax = 5°/Myr)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "lcurve.png"), dpi=130, bbox_inches="tight")
print("wrote lcurve.png")
