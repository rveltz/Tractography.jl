import ForwardDiff

"""
$(SIGNATURES)

Get the ODF vector length corresponding to a given lmax.
This length is the result of

```
n = 0
for l = 0:2:lmax, m = -l:l
    n += 1
end
n
# gives sequence 1  6  15  28  45  66  91  120
```
"""
get_odf_length(lmax) = div(lmax^2,2) + div(3*lmax,2) + 1

"""
$(SIGNATURES)

Get the lmax from the ODF length, ie `size(data, 4)`. The is the inverse mapping of `get_odf_length`.
"""
get_lmax_from_odf_length(n) = Int(-3/2 + sqrt(1 + 8*n)/2)

"""
$(SIGNATURES)

Evaluate the θ derivative of real orthonormal spherical harmonics.

Based on `FastTransforms.sphevaluate` and `ForwardDiff.jl`
"""
function ∂θro_sh(θ, ϕ, l, m, outer_f = identity)
    ForwardDiff.derivative(z -> outer_f(ro_sh(z, ϕ, l, m)), θ)
end

"""
$(SIGNATURES)

Evaluate the ϕ derivative of real orthonormal spherical harmonics.

Based on `FastTransforms.sphevaluate` and `ForwardDiff.jl`
"""
function ∂ϕro_sh(θ, ϕ, l, m, outer_f = identity)
    ForwardDiff.derivative(z -> outer_f(ro_sh(θ, z, l, m)), ϕ)
end

"""
$(SIGNATURES)

Evaluate the real orthonormal spherical harmonics.

Recall that ∫ ODF = c₀₀ √4π

Based on `FastTransforms.sphevaluate`. You can also check the [url](https://juliaapproximation.github.io/FastTransforms.jl/dev/#FastTransforms.sphevaluate)
"""
function ro_sh(θ, ϕ, l, m)
    FastTransforms.sphevaluate(θ, ϕ, l, m)
end

"""
$(SIGNATURES)

Evaluate the real orthonormal spherical harmonics Yₗₘ(θ, ϕ) for `(θ, Φ) ∈ angles` for l ∈ [0, lₘₐₓ] and m ∈ [-l, l]. 

It returns a 2d array `Yₗₘ[(θ, ϕ), index]`.

It is mainly used to cache the harmonics Yₗₘ for later use (evaluation of ODF).

# Arguments
- `angles` a vector of tuples
- `lmax::Int` maximum l for harmonics.

> The storage convention is explained in https://mrtrix.readthedocs.io/en/dev/concepts/spherical_harmonics.html#storage-conventions.
"""
function get_vector_of_sh(angles::AbstractVector{Tuple{𝒯, 𝒯}}, lmax, der::Int = 0; outer_f = identity) where 𝒯
    odf_length = get_odf_length(lmax)
    Yₗₘ = zeros(𝒯, length(angles), odf_length)
    for (i, (θ, ϕ) ) in pairs(angles)
        n = 1
        for l = 0:2:lmax, m = -l:l
            if der == 0
                y = ro_sh(θ, ϕ, l, m)
            elseif der == 1
                y = ∂θro_sh(θ, ϕ, l, m, outer_f)
            else
                y = ∂ϕro_sh(θ, ϕ, l, m, outer_f)
            end
            Yₗₘ[i, n] = y * (-1)^m
            n += 1
        end
    end
    Yₗₘ
end