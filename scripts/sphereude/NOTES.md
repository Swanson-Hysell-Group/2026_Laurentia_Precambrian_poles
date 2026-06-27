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

## Current choice (provisional)

`λ1 = 3e5`, `ωmax = 2.5°/Myr`, `niter = 2000` (ADAM + L-BFGS). From a quick
three-point scan at `ωmax = 2.5`:

| λ1   | total arc | peak rate  | character                                  |
|------|-----------|-----------|--------------------------------------------|
| 3e4  | ~439°     | 1.75°/Myr | wiggly; chases poles, loops in Keweenawan   |
| 3e5  | ~377°     | 1.27°/Myr | balanced; threads knot, bends to old end    |
| 3e6  | ~287°     | 0.73°/Myr | very schematic; shortcuts old end, slow     |

`3e5` keeps a peak rate above 1°/Myr (consistent with rapid Laurentia motion) and
recovers the sparse old end, without the coiling of `3e4`.

## Next steps / to revisit

1. **Pick `λ1` by L-curve elbow.** Sweep `λ1` (e.g. 1e4, 3e4, 1e5, 2e5, 3e5, 5e5,
   1e6, 3e6, 1e7) at `ωmax = 2.5`, and for each record the empirical misfit
   `Σκ(1−cosθ)` against a roughness proxy `Σ‖d²x‖²`. Plot log–log and take the
   elbow (Hansen, 2001; Gallo et al., 2022 elbow strategy). `analyze_lcurve.py`
   in this folder computes the misfit, roughness, peak rate, and arc length for a
   set of `lcurve_<λ1>.csv` fits and plots the L-curve; it expects the sweep
   outputs from running `fit_apwp_spline.jl` with `SU_LAMBDA1`/`SU_OUT` set per λ1.
2. **Likely range is [3e4, 3e5].** We may prefer something slightly looser than 3e5
   to let the path honor more structure — to be chosen against the L-curve rather
   than by eye.
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
