# SphereUDE APWP spline fit

This directory is a **self-contained Julia project** that fits a spherical
apparent-polar-wander path to the Laurentia poles using
[SphereUDE.jl](https://github.com/ODINN-SciML/SphereUDE.jl) (Sapienza et al.,
2025, *JGR Machine Learning and Computation*, doi:10.1029/2025JH000626).

The path is the solution of `dx/dt = L(t) × x(t)`, where the time-dependent Euler
vector `L(t)` is learned by a neural network and regularized to be slow and smooth
(order `k = 1`, power `p = 2`), following the paper's real-data recipe.

## Two ways to reproduce

**1. Without Julia (figures only).** The fitted path is committed to
`data/nordic_summaries/apwp_sphereude_path.csv`. `scripts/build_apwp_figure.py`
reads that CSV directly and overlays the path, so the figures reproduce from the
Python environment alone — no Julia required. If the CSV is absent the figures
are simply drawn without the path.

**2. Re-running the fit (requires Julia).** The Julia dependencies are pinned
**locally to this directory** in `Project.toml` / `Manifest.toml`; they are kept
out of the project's Python environment on purpose.

```
# one-time: instantiate the pinned environment
julia +lts --project=scripts/sphereude -e 'using Pkg; Pkg.instantiate()'

# export the fit input (rotated, filtered poles), then fit, then re-plot
python scripts/build_apwp_figure.py            # writes data/.../apwp_fit_input.csv
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
| `SU_LAMBDA1`  | `1e6`   | penalty on `‖dL/dt‖²` — higher → smoother path       |
| `SU_LAMBDA0`  | `1e0`   | penalty on `‖L‖²` — higher → slower path             |
| `SU_OMEGADEG` | `2.5`   | cap on angular velocity, deg/Myr (soft rate bound)  |
| `SU_NITER`    | `2000`  | ADAM and L-BFGS iterations each                     |
| `SU_OUT`      | path CSV| output CSV location                                 |

`ωmax` is set as a *loose physical* bound (2.5 deg/Myr): it permits the rapid
(>1 deg/Myr) Keweenawan motion recorded by the poles while preventing runaway
excursions. The path smoothness is controlled by `λ1`, chosen by L-curve analysis.

`build_apwp_figure.py` excludes two groups of poles from the fit input (both are
still plotted): the ca. 1382 Ma Greenland poles (Midsommersø, Victoria Fjord,
Zig-Zag Dal), which conflict with the rest of the path, and any pole whose age
half-range exceeds 50 Myr, which drops the loosely-dated Scotland poles (Torridon,
Stoer) so they do not pull on the path.
