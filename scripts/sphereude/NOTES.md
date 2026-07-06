# Notes: fitting the Laurentia APWP spline with SphereUDE

Working notes on how the path in `data/nordic_summaries/apwp_sphereude_path.csv`
was produced, the modeling choices behind it, and what to revisit. See
`README.md` for how to run things; this file is the *why* and the *to-do*.

## Method

We fit a single spherical apparent-polar-wander path to the Laurentia poles with
SphereUDE.jl (Sapienza et al., 2025, doi:10.1029/2025JH000626). The path solves
`dx/dt = L(t) × x(t)`, where the time-dependent Euler vector `L(t)` (rotation axis
+ angular rate) is parametrized by a neural network and learned by minimizing a
von Mises–Fisher data misfit plus regularization. Poles are weighted by their
Fisher concentration `κ ≈ (140°/A95)²` (McElhinny & McFadden, 1999).

## Data selection (done in `build_apwp_figure.py::export_fit_input`)

The fit input is the rotated, age-restricted pole set (717–1779 Ma) with two
exclusions; both excluded groups are still *plotted*, just not fit:

1. **ca. 1382 Ma Greenland poles** (Midsommersø, Victoria Fjord, Zig-Zag Dal) —
   conflict with the rest of the path.
2. **Poles with age half-range > 50 Myr** — loosely-dated poles should not pull on
   a path that is parametrized by age. At 50 Myr this drops exactly the two
   Scotland poles (Torridon ±145 Myr, Stoer ±70 Myr) and nothing else (next
   loosest is Adirondack at ±23 Myr).

Net: 53 poles enter the fit.

## The two regularization knobs

- **`ωmax` (rate cap, currently 2.5°/Myr).** A hard bound on `‖L(t)‖`, the angular
  velocity. Set as a *loose physical* regularizer: it permits the rapid
  (>1°/Myr) Keweenawan motion that the poles record, while preventing runaway
  excursions. It is deliberately NOT tight — capping near 1°/Myr would manufacture
  slow rates by construction and contradict the rapid-motion interpretation.
- **`λ1` (smoothness, `SU_LAMBDA1`).** Coefficient on `∫‖dL/dt‖² dt`: penalizes how
  fast `L(t)` itself changes, i.e. how abruptly the direction/rate of motion turns.
  This is the main bias–variance lever and is applied globally over the age span.
  - Too low → overfit: the path coils and chases individual poles.
  - Too high → underfit: the path takes smooth shortcuts and misses isolated
    poles where data are sparse (e.g. the ~1590/1635 Ma old end), and it
    suppresses the real rapid intervals.

A rapid but *sustained* rotation has large `‖L‖` yet small `‖dL/dt‖`, so a generous
`ωmax` plus a moderate `λ1` can capture fast Keweenawan motion while staying
schematic — the rate cap and the smoothness penalty are doing different jobs.

## Current choice

`λ1 = 2e5`, `ωmax = 2.5°/Myr`, `niter = 2000` (ADAM + L-BFGS) is the committed
default and the path in `apwp_sphereude_path.csv` / the main Lambert figure. From
the sweep (`compare_lambda.py`, `lambda_sweep/`):

| λ1   | misfit | roughness | peak rate  | note                                  |
|------|--------|-----------|-----------|----------------------------------------|
| 3e4  | 297    | 2.41e-4   | 1.75°/Myr | wiggly; chases poles, loops             |
| 5e4  | 340    | 1.85e-4   | 1.69°/Myr |                                         |
| 7e4  | 378    | 1.53e-4   | 1.62°/Myr |                                         |
| 1e5  | 427    | 1.28e-4   | 1.52°/Myr | **L-curve elbow** (most faithful balance)|
| 2e5  | 526    | 1.06e-4   | 1.37°/Myr | **committed** (minimum roughness)        |
| 3e5  | 582    | 1.18e-4   | 1.27°/Myr | roughness turns back up (over-smoothed)  |

We chose `2e5` because it **minimizes path roughness** — the smoothest, most
schematic path — which suits presenting the path as schematic rather than
definitive, while still keeping a peak rate above 1°/Myr and reaching the sparse
old end. `1e5` is the L-curve elbow and the more data-faithful alternative (it
reaches the ca. 1140 Ma apex and cuts more directly to the ca. 990 Ma poles); it
is the natural choice if we later want the more quantitative path.

Note the committed CSV was produced at `niter = 1500` (from the sweep); rerunning
the default (`niter = 2000`) gives a near-identical but not bit-identical path.

## Next steps / to revisit

1. **Pick `λ1` by L-curve elbow.** A sweep at `ωmax = 2.5` over
   {3e4, 5e4, 7e4, 1e5, 2e5, 3e5} is saved in `lambda_sweep/elbow_<λ1>.csv`.
   `compare_lambda.py` reads those, prints misfit / roughness / peak rate / arc
   length, and writes two figures to `_static/`:
   `Laurentia_apwp_lambda_lcurve.png` (roughness vs. misfit) and
   `Laurentia_apwp_lambda_comparison.png` (one map panel per λ1). The L-curve
   **elbow sits at λ1 ≈ 1e5** (misfit 427, roughness 1.28e-4, peak 1.52°/Myr):
   below it roughness rises steeply for little misfit gain; above it (2e5, 3e5)
   misfit climbs with little further smoothing.
   To regenerate the sweep (requires Julia), fit each λ1 with `SU_LAMBDA1`/`SU_OUT`
   pointed at `lambda_sweep/` (see the header of `compare_lambda.py`).
2. **Committed = 2e5 (min roughness); revisit 1e5 (elbow) if a more faithful path
   is wanted.** 2e5 is the smoothest/most schematic; 1e5 is the L-curve elbow that
   reaches the ca. 1140 Ma apex and cuts more directly to the ca. 990 Ma poles.
   Swapping is a one-line default change plus a rerun (see `Current choice`).
3. **Report the peak-rate-vs-age profile** of the chosen fit so the rapid intervals
   (Keweenawan) are quantified as a result, not assumed.
4. **Sensitivity checks:** vary the random seed and `niter` to confirm the path is
   stable (the paper notes inversions can depend on NN architecture/seed); consider
   averaging over runs or reporting spread.
5. **Held in reserve — split fit.** If a single global `λ1` cannot be both smooth
   through the Keweenawan and faithful elsewhere, fit separate splines before and
   after the Keweenawan rapid interval and join them. This is a fallback, not the
   first option.
6. **Uncertainty quantification.** SphereUDE's roadmap (and Schmid et al., 2025)
   points to Bayesian/ensemble uncertainty on the path; revisit if we want a
   confidence envelope rather than a single schematic curve.
