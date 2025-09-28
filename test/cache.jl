using Test, LinearAlgebra, Accessors
using Tractography
const TG = Tractography

model = TMC(Δt = 0.125f0,
            odfdata = ODFData("../examples/fod-FC.nii.gz"),
            C = Cone(45),
            )

alg = Deterministic()
cache_c = TG._init(model, (CSD()); n_sphere = 200)

@test all(isfinite, cache_c.odf)

# test norm of directions
@test all(x -> abs(norm(x) - 1) < 1e-6, eachrow(cache_c.directions))

# test spherical coordinates
@test all(x -> (0 <= x[1] <= Float32(pi)),  cache_c.angles)
@test all(x -> (0 <= x[2] <= 2pi), cache_c.angles)