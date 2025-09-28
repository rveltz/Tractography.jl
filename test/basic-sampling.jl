using Test, Tractography, LinearAlgebra
Tractography.spherical_to_euclidean(0,0)
@test all(Tractography.euclidean_to_spherical(Tractography.spherical_to_euclidean(0.1, -0.01)...) .≈ [0.1, -0.01])
u0 = normalize(rand(3))
@test all([Tractography.spherical_to_euclidean(Tractography.euclidean_to_spherical(u0...)...)...] .≈ u0)


p0 = normalize(rand(3))
v0 = normalize(rand(3)); v0 .-= dot(p0,v0) .* p0
@test dot(v0, p0) ≈ 0 atol = 1e-14
