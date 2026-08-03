# Fit a spherical apparent-polar-wander path to the Laurentia poles with
# SphereUDE.jl (Sapienza et al., 2025, JGR Machine Learning and Computation,
# doi:10.1029/2025JH000626). The path is the solution of dx/dt = L(t) x x(t)
# with the time-dependent Euler vector L(t) learned by a neural network and
# regularized to be slow and smooth (order k = 1, power p = 2), following the
# paper's real-data (Gondwana) recipe.
#
# Inputs are written by build_apwp_figure.py (age, plat, plon, a95; rotated into
# Laurentia coordinates and filtered -- Svalbard and the ca. 1382 Ma Greenland
# poles are already excluded). Run it first to (re)generate them.
#
# Two pole sets are fit and kept separately, differing only in how the
# sedimentary poles enter. SphereUDE weights each pole by a single Fisher
# concentration kappa ~ (140/A95)^2, so a Kent ellipse cannot be consumed
# directly; it enters the corrected track only through its equal-area-equivalent
# scalar A95 = sqrt(zeta95 * eta95), discarding any directionality.
#
#   apwp_fit_input_corrected.csv   -> apwp_sphereude_path_corrected.csv
#       sedimentary poles at their inclination-shallowing-corrected (Kent mean)
#       positions. This is the default here, and the track the figures overlay.
#   apwp_fit_input_uncorrected.csv -> apwp_sphereude_path_uncorrected.csv
#       the same poles as measured, without the correction.
#
# Defaults fit the corrected track:
#
#   julia +lts --project=scripts/sphereude scripts/sphereude/fit_apwp_spline.jl
#
# The uncorrected track, via SU_IN / SU_OUT:
#
#   SU_IN=data/nordic_summaries/apwp_fit_input_uncorrected.csv \
#   SU_OUT=data/nordic_summaries/apwp_sphereude_path_uncorrected.csv \
#     julia +lts --project=scripts/sphereude scripts/sphereude/fit_apwp_spline.jl
#
# Use the 1.10 LTS: the SciML/Zygote stack segfaults during gradient JIT on Julia
# 1.12, and Manifest.toml is resolved for 1.10 (see README.md).

using SphereUDE
using SciMLSensitivity
using Optimization, OptimizationOptimisers, OptimizationOptimJL
using CSV, DataFrames
using LinearAlgebra, Random

const HERE = @__DIR__
const ROOT = dirname(dirname(HERE))
const IN_CSV = get(ENV, "SU_IN",
    joinpath(ROOT, "data", "nordic_summaries", "apwp_fit_input_corrected.csv"))
const OUT_CSV = get(ENV, "SU_OUT",
    joinpath(ROOT, "data", "nordic_summaries", "apwp_sphereude_path_corrected.csv"))

# Hyperparameters (overridable via environment for tuning sweeps):
#   SU_LAMBDA1  smoothness penalty on dL/dt (higher -> more schematic)
#   SU_LAMBDA0  magnitude penalty on |L|    (higher -> slower path)
#   SU_OMEGADEG cap on angular velocity, deg/Myr
#   SU_NITER    iterations for each of ADAM and LBFGS
# λ1 = 2e5 (see NOTES.md): the value that minimizes path roughness in the sweep,
# giving the smoothest, most schematic path -- appropriate since the path is
# presented as schematic rather than definitive. Peak rate ~1.4 deg/Myr (>1,
# consistent with rapid Laurentia motion) under the 2.5 deg/Myr cap. The L-curve
# elbow (more faithful balance) is nearer 1e5; compare_lambda.py shows the sweep.
const LAMBDA1 = parse(Float64, get(ENV, "SU_LAMBDA1", "2.0e5"))
const LAMBDA0 = parse(Float64, get(ENV, "SU_LAMBDA0", "1.0e0"))
const OMEGADEG = parse(Float64, get(ENV, "SU_OMEGADEG", "2.5"))
const NITER = parse(Int, get(ENV, "SU_NITER", "2000"))

rng = Random.default_rng()
Random.seed!(rng, 1234)

df = CSV.read(IN_CSV, DataFrame)
sort!(df, :age)
println("Fitting SphereUDE path to $(nrow(df)) poles, $(minimum(df.age))-$(maximum(df.age)) Ma")

# pole lat/lon (degrees) -> 3 x N matrix of unit vectors
X = sph2cart(permutedims(hcat(df.plat, df.plon)); radians = false)

# A95 -> Fisher concentration kappa (McElhinny & McFadden, 1999): A95 ~ 140/sqrt(kappa)
kappas = (140.0 ./ Float64.(df.a95)) .^ 2

times = Float64.(df.age)
tmin, tmax = minimum(times), maximum(times)

data = SphereData(times = times, directions = X, kappas = kappas, L = nothing)

# slow + smooth regularization on L(t): penalize its magnitude (order 0) and its
# first derivative (order 1, power 2), as in the paper's Gondwana fit.
regs = [
    Regularization(order = 1, power = 2.0, λ = LAMBDA1, diff_mode = FiniteDiff(1.0e-6)),
    Regularization(order = 0, power = 2.0, λ = LAMBDA0, diff_mode = nothing),
]
println("Hyperparameters: λ1=$LAMBDA1  λ0=$LAMBDA0  ωmax=$OMEGADEG deg/Myr  niter=$NITER")

params = SphereParameters(
    tmin = tmin,
    tmax = tmax,
    reg = regs,
    train_initial_condition = true,   # Precambrian path: u0 is free, not the spin axis
    multiple_shooting = false,
    pretrain = false,
    u0 = [0.0, 0.0, 1.0],
    ωmax = OMEGADEG * π / 180.0,       # cap angular velocity (deg/Myr -> rad/Myr)
    reltol = 1.0e-7,
    abstol = 1.0e-7,
    niter_ADAM = NITER,
    niter_LBFGS = NITER,
    verbose_step = 500,
    sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP(true)),
)

results = train(data, params, rng, nothing, nothing)
println("Final loss: ", results.losses[end])

# results.fit_directions is 3 x 1000 unit vectors at results.fit_times; convert to lat/lon
fit_sph = cart2sph(results.fit_directions; radians = false)   # 2 x 1000 -> [lat; lon]
lat = fit_sph[1, :]
lon = mod.(fit_sph[2, :], 360.0)                              # 0-360 to match pole longitudes

out = DataFrame(age = results.fit_times, lat = lat, lon = lon)
CSV.write(OUT_CSV, out)
println("Wrote ", OUT_CSV, " (", nrow(out), " path samples)")
