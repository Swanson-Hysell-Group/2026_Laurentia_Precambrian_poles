# Notes: fitting the Laurentia APWP spline with SphereUDE

Working notes on how the paths in
`data/nordic_summaries/apwp_sphereude_path_{corrected,uncorrected}.csv` were
produced, the modeling choices behind it, and what to revisit. See
`README.md` for how to run things; this file is the *why* and the *to-do*.

## Two tracks: corrected and uncorrected poles

`build_apwp_figure.py` applies `apply_kent_poles` inside `load_path`, so the
sedimentary units are substituted with their inclination-shallowing-corrected
(Kent mean) poles. Two pole sets are therefore fit and kept separately:

| track | input | output | sedimentary poles |
|-------|-------|--------|-------------------|
| `corrected` | `apwp_fit_input_corrected.csv` | `apwp_sphereude_path_corrected.csv` | Kent-corrected |
| `uncorrected` | `apwp_fit_input_uncorrected.csv` | `apwp_sphereude_path_uncorrected.csv` | as measured |

The corrected set is the pipeline default (`FIGURE_TRACK = "corrected"` in
`build_apwp_figure.py`), and is what the map figures and
`build_reconstruction_figure.py` use, so the overlaid path and the plotted poles
come from the same pole positions. Both inputs are regenerated on every
`build_apwp_figure.py` run; neither is a hand-kept snapshot.

`load_spline_path` reads each track only from its own file, with no fallback to
the pre-track names `apwp_fit_input.csv` / `apwp_sphereude_path.csv`. Those names
are not self-describing, and a fallback silently mislabels whichever pole set the
file happens to hold -- which is exactly what happened on 2026-07-31, when a
corrected fit written to the legacy name was served as the uncorrected track.
The two legacy files are superseded and referenced by nothing.

A note on the naming: the tracks are named by whether the correction has been
applied, rather than "flattened"/"unflattened". Compaction shallows (flattens)
sedimentary inclinations and unflattening *is* the correction, so "unflattened"
denotes the *corrected* set — a pair that inverts easily in reading. Note that
`pole_tools.make_nordic_summary`'s `pole_mean_unflattened` argument and
`unflattened_pole` in the Lower Freda notebook still use that older vocabulary
for the corrected pole.

Pass `SU_IN` / `SU_OUT` to fit either track; see the header of
`fit_apwp_spline.jl` for the exact invocation.

**Why the correction cannot be carried in full.** SphereUDE weights each pole by
a single Fisher concentration `κ ≈ (140°/A95)²`, so it has no way to represent a
Kent ellipse's two concentration parameters. The ellipse enters only through the
equal-area-equivalent circular radius `A95 = sqrt(ζ95 · η95)` that
`apply_kent_poles` computes — an isotropic simplification of an anisotropic
uncertainty. Any directionality in the Kent confidence regions is discarded by
the fit, which is worth stating whenever the corrected path is shown.

The correction is not cosmetic: it moves 11 of the 54 poles, the largest being
Torridon at 13.8°, and inflates their A95 (mean 4.4° → 6.1° over the moved
poles), so those poles also pull less on the corrected fit.

### Both tracks fit 2026-08-01

Both at the committed defaults (λ1 = 2e5, λ0 = 1, ωmax = 2.5°/Myr,
`niter = 2000`), against inputs that include the revised ca. 1108 Ma Nipigon pole
and the final Kent means, so the two are directly comparable.

| | uncorrected | corrected |
|---|---|---|
| final loss | 32.88 | 30.85 |
| empirical term | 23.86 | 22.39 |
| order-1 reg | 8.90 | 8.36 |
| peak rate | 1.40°/Myr at 1098 Ma | 1.32°/Myr at 1096 Ma |
| arc length | 401.9° | 367.3° |
| mean angular misfit to own poles | 5.58° | 5.97° |
| max misfit | 14.1° at 1109 Ma | 15.0° at 1109 Ma |

The two paths separate by 2.7° on average (median 1.6°), with a maximum of 13.9°
at ca. 1380 Ma — near the Belt Supergroup sedimentary poles (McNamara,
Pilcher/Garnet Range) and the ca. 975 Ma Torridon pole, which the correction
moves most.

Two points worth care when reporting these:

1. **The corrected fit's lower loss is not evidence of a better fit.** Its
   *angular* misfit to its own poles is slightly worse (5.97° vs 5.58°); the loss
   falls because the Kent A95s are larger, so `κ` is smaller and each pole is
   down-weighted. Loss is comparable across λ1 at fixed data, not across data
   sets with different uncertainties.
2. **The peak rate is essentially unchanged** (1.32 vs 1.40°/Myr, both at ca.
   1096–1098 Ma). The Keweenawan rapid-motion interval is carried by igneous
   poles, which the correction does not touch, so this is the expected result and
   is a useful check that the correction did not leak into that part of the path.
   The corrected path is ~35° shorter overall, i.e. less sinuous.

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
   a path that is parametrized by age. Stoer (±70 Myr) is now the *only* pole in
   the compilation this drops. The next loosest admitted is Torridon at exactly
   ±50.0 Myr (see below), and after that Adirondack at ±23 Myr — nothing in the
   file falls between 23 and 70, so the threshold separates no natural
   population; it passes straight through the single pole sitting at its value.

Net: 54 poles enter the fit (as of the 2026-07-31 rerun).

**Torridon now enters the fit.** The Torridon age update narrowed its range to
950–1050 Ma, i.e. a half-range of exactly 50.0 Myr. The filter is strictly
`> FIT_MAX_AGE_HALF_UNC`, so Torridon passes by a knife-edge and is now fit as a
ca. 975 Ma pole (rotated: −21.8°N / 188.9°E, A95 7.1 uncorrected; −33.0°N /
179.8°E, A95 10.7 corrected — the single largest pole shift the Kent correction
makes). It is the oldest-end
Scotland constraint on the ca. 990 Ma part of the path, and it is worth being
aware that a 1 Myr widening of its age range would silently drop it again.
Chengwatana (ca. 1096 Ma) left the fit because it was removed upstream from
`nordic_summaries_combined.csv`, not by any filter here.

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
default for both tracks and the path in the main Lambert figure. From
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

**Refit 2026-07-31, uncorrected track** (`apwp_sphereude_path_uncorrected.csv`)
against the updated pole set (54 poles, 719–1779 Ma) at the committed defaults
and `niter = 2000` (2000 ADAM + 2000 L-BFGS). Final loss 32.39
(empirical 23.87, order-0 reg 0.114, order-1 reg 8.41). The path is stable against
the constraint changes: mean displacement from the previous committed path is
1.8°, median 1.0°, max 8.0° near 1201 Ma. Peak rate is 1.37°/Myr at ca. 1097 Ma —
unchanged from the λ1 = 2e5 sweep entry below — and arc length is 401°.

The sweep table below was computed on the *previous* pole set and its misfit /
roughness columns are therefore stale as absolute numbers; it is retained because
the *ordering* and the location of the elbow are what the choice rests on. Rerun
`compare_lambda.py` over a fresh sweep before quoting those values.

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
