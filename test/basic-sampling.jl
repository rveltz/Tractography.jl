using Test, Tractography, LinearAlgebra, Accessors

Tractography.PlottingSH()
Tractography.PreComputeAllODF()
Tractography.DirectSH()

Tractography.softplus(0,1)
Tractography.∂softplus(0,1)

Tractography.spherical_to_euclidean(0,0)
@test all(Tractography.euclidean_to_spherical(Tractography.spherical_to_euclidean(0.1, -0.01)...) .≈ [0.1, -0.01])
u0 = normalize(rand(3))
@test all([Tractography.spherical_to_euclidean(Tractography.euclidean_to_spherical(u0...)...)...] .≈ u0)


p0 = normalize(rand(3))
v0 = normalize(rand(3)); v0 .-= dot(p0,v0) .* p0
@test dot(v0, p0) ≈ 0 atol = 1e-14
# u = Tractography.Exp𝕊²(p0, v0, 0.2)
# @test norm(u) ≈ 1

model = Tractography.TMC(Δt = 0.125f0,
            odfdata = Tractography.ODFData((@__DIR__) * "/../examples/fod-FC.nii.gz"),
            cone = Tractography.Cone(15),
            proba_min = 0.005f0,
            )

model_diffusion = Tractography.TMC(Δt = 0.001f0,
            odfdata = Tractography.ODFData((@__DIR__) * "/../examples/fod-FC.nii.gz"),
            proba_min = 0.0f0,
            evaluation_algo = Tractography.DirectSH(),
            )

Tractography._apply_mask!(model, ones(64, 64, 3))

show(stdout, model)

Tractography.sample(model, Tractography.Deterministic(), rand(Float32, 6, 2); nt = 10, maxodf_start = true, reverse_direction = true);
Tractography.sample(model, Tractography.Connectivity(Tractography.Deterministic()), rand(Float32, 6, 2); nt = 10, maxodf_start = true);
Tractography.sample(model, Tractography.Probabilistic(), rand(Float32, 6, 2); nt = 10);

Tractography.sample(model_diffusion, Tractography.Diffusion(), rand(Float32, 6, 2); nt = 10);
Tractography.sample(model, Tractography.Diffusion(), rand(Float32, 6, 2); nt = 10);

########################
# cache
Nmc = 10
seeds = rand(Float32, 6, Nmc)
streamlines = zeros(Float32, 6, 20, Nmc)
tract_length = zeros(UInt32, Nmc)
alg = Probabilistic()
cache = Tractography.init(model, alg)
show(stdout, cache)
cache = Tractography._init(model, alg)
show(stdout, cache)

########################
# cache diffusion
show(Tractography.Diffusion())
cache = Tractography.init(model_diffusion, Tractography.Diffusion())
cache = Tractography.init((@set model_diffusion.evaluation_algo = Tractography.PreComputeAllODF()), Tractography.Diffusion())
