# SphereUDE APWP spline fit

This directory is a **self-contained Julia project** that fits a spherical
apparent-polar-wander path to the Laurentia poles using
[SphereUDE.jl](https://github.com/ODINN-SciML/SphereUDE.jl) (Sapienza et al.,
2025, *JGR Machine Learning and Computation*, doi:10.1029/2025JH000626).

The path is the solution of `dx/dt = L(t) × x(t)`, where the time-dependent Euler
vector `L(t)` is learned by a neural network and regularized to be slow and smooth
(order `k = 1`, power `p = 2`), following the paper's real-data recipe.

## Two ways to reproduce

**1. Without Julia (figures only).** The fitted paths are committed to
`data/nordic_summaries/apwp_sphereude_path_corrected.csv` and
`..._uncorrected.csv`. `scripts/build_apwp_figure.py` reads them directly and
overlays the `corrected` one (`FIGURE_TRACK`), so the figures reproduce from the
Python environment alone — no Julia required. If the CSV is absent the figures
are simply drawn without the path.

**Two tracks.** The paths differ only in how the sedimentary poles enter:
`corrected` uses the inclination-shallowing-corrected (Kent mean) positions,
`uncorrected` takes them as measured. The figures plot Kent-corrected poles, so
they overlay the `corrected` path. See `NOTES.md` for the comparison and for why
the Kent ellipse can only enter as its equal-area circular `A95`.

**2. Re-running the fit (requires Julia).** The Julia dependencies are pinned
**locally to this directory** in `Project.toml` / `Manifest.toml`; they are kept
out of the project's Python environment on purpose.

```
# one-time: instantiate the pinned environment
julia +lts --project=scripts/sphereude -e 'using Pkg; Pkg.instantiate()'

# export both fit inputs (rotated, filtered poles), then fit, then re-plot
python scripts/build_apwp_figure.py     # writes apwp_fit_input_{corrected,uncorrected}.csv

# corrected track (the default, and what the figures overlay)
julia +lts --project=scripts/sphereude scripts/sphereude/fit_apwp_spline.jl

# uncorrected track
SU_IN=data/nordic_summaries/apwp_fit_input_uncorrected.csv \
SU_OUT=data/nordic_summaries/apwp_sphereude_path_uncorrected.csv \
  julia +lts --project=scripts/sphereude scripts/sphereude/fit_apwp_spline.jl

python scripts/build_apwp_figure.py            # overlays the refit path
```

## Julia version

Use the **1.10 LTS** (`julia +lts`, via [juliaup](https://github.com/JuliaLang/juliaup)).
The SciML/Zygote stack SphereUDE depends on segfaults during gradient JIT on
Julia 1.12; `Manifest.toml` is resolved for 1.10.

## Hyperparameters

Set as defaults in `fit_apwp_spline.jl`, overridable via environment variables for
tuning sweeps:

| env var       | default | meaning                                             |
|---------------|---------|-----------------------------------------------------|
| `SU_LAMBDA1`  | `2e5`   | penalty on `‖dL/dt‖²` — higher → smoother path       |
| `SU_LAMBDA0`  | `1e0`   | penalty on `‖L‖²` — higher → slower path             |
| `SU_OMEGADEG` | `2.5`   | cap on angular velocity, deg/Myr (soft rate bound)  |
| `SU_NITER`    | `2000`  | ADAM and L-BFGS iterations each                     |
| `SU_IN`       | corrected| input CSV location                                 |
| `SU_OUT`      | corrected| output CSV location                                |

`ωmax` is set as a *loose physical* bound (2.5 deg/Myr): it permits the rapid
(>1 deg/Myr) Keweenawan motion recorded by the poles while preventing runaway
excursions. The path smoothness is controlled by `λ1`, chosen by L-curve analysis.

`build_apwp_figure.py` excludes two groups of poles from the fit input (both are
still plotted): the ca. 1382 Ma Greenland poles (Midsommersø, Victoria Fjord,
Zig-Zag Dal), which conflict with the rest of the path, and any pole whose age
half-range exceeds 50 Myr, which drops the loosely-dated Scotland poles (Torridon,
Stoer) so they do not pull on the path.
