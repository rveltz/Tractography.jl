"""
$(TYPEDEF)

Cache for computing streamlines in `sample`.

# Arguments (with default values):

$(TYPEDFIELDS)
"""
@with_kw_noshow struct TMCCache{𝒯y, 𝒯a, 𝒯dir, 𝒯d, 𝒯c, 𝒯all, 𝒯alld1, 𝒯alld2} <: AbstractCache
    "Matrix for holding the real orthonormal spherical harmonics sampled on a grid."
    Yₗₘ::𝒯y = nothing
    "Matrix for holding the θ derivative of orthonormal spherical harmonics sampled on a grid."
    ∂θYₗₘ::𝒯y = Yₗₘ
    "Matrix for holding the ϕ derivative of orthonormal spherical harmonics sampled on a grid."
    ∂ϕYₗₘ::𝒯y = Yₗₘ
    "List of angles (θ, ϕ) for sampling the sphere."
    angles::𝒯a = nothing
    "List of directions ∈ 𝕊² for sampling the sphere."
    directions::𝒯dir = nothing
    "lmax used for real orthonormal spherical harmonics."
    lmax::Int
    "Number of points in the directions."
    n_sphere::Int
    "Measure to compute integrals of probabilities."
    dΩ::𝒯d
    "Buffer to hold the cone sample on the grid."
    cone::𝒯c = nothing
    "Matrix of all ODF values."
    odf::𝒯all = nothing
    "Matrix of all ∂θ ODF values."
    ∂θodf::𝒯alld1 = nothing
    "Matrix of all ∂ϕ ODF values."
    ∂ϕodf::𝒯alld2 = nothing
end
@inline get_angles(cache, ii) = cache.angles[ii]

function Base.show(io::IO, cache::TMCCache{𝒯y, 𝒯a, 𝒯dir, 𝒯d, 𝒯c, 𝒯all, 𝒯alld1, 𝒯alld2}) where {𝒯y, 𝒯a, 𝒯dir, 𝒯d, 𝒯c, 𝒯all, 𝒯alld1, 𝒯alld2}
    printstyled(io, "TMCCache (subset)"; bold = true, color = :cyan)
    printstyled(io, "\n ├─ size : $(round(Base.summarysize(cache)/1024^3, digits=3)) GiB"; bold = true)
    println(io,   "\n ├─ dΩ::$Td  : ", cache.dΩ)
    println(io,   " ├─ Yₗₘ  : ", 𝒯y)
    println(io,   " ├─ ∂θodf : ", 𝒯alld1)
    if ~isnothing(cache.∂θodf)
        println(io, "      └─── size = ", size(cache.∂θodf))
    end
    println(io,   " ├─ ∂ϕodf : ", 𝒯alld2)
    if ~isnothing(cache.∂ϕodf)
        println(io, "      └─── size = ", size(cache.∂ϕodf))
    end
    println(io,   " └─ odf   : ", 𝒯all)
    if ~isnothing(cache.odf)
        println(io, "      └─── size = ", size(cache.odf))
    end
end

"""
$(TYPEDEF)

Cache specific to threaded or GPU computations.

# Fields
$(TYPEDFIELDS)
"""
struct ThreadedCache{𝒯a, 𝒯c, 𝒯ai, 𝒯ang, 𝒯, 𝒯s} <: AbstractCache
    "Matrix of all ODF values."
    odf::𝒯a
    "Matrix of all ∂θ ODF values."
    ∂θodf::Union{Nothing, 𝒯a}
    "Matrix of all ∂ϕ ODF values."
    ∂ϕodf::Union{Nothing, 𝒯a}
    "Buffer to hold the cone sample on the grid."
    cone::𝒯c
    "List of directions ∈ 𝕊² for sampling the sphere."
    directions::𝒯ai
    "List of angles (θ, ϕ) for sampling the sphere."
    angles::𝒯ang
    "Integral of fodf (after correction with mollifier)."
    ∫odf::𝒯s
    "Measure to compute integrals of probabilities."
    dΩ::𝒯
end

function Base.show(io::IO, cache::ThreadedCache{𝒯a}) where {𝒯a}
    printstyled(io, "ThreadedCache (subset)"; bold = true, color = :cyan)
    printstyled(io, "\n ├─ size : $(round(Base.summarysize(cache)/1024^3, digits=3)) GiB\n"; bold = true)
    println(io,   " ├─ ∂θodf : ", typeof(cache.∂θodf))
    if ~isnothing(cache.∂θodf)
        println(io, "      └─── size = ", size(cache.∂θodf))
    end
    println(io,   " ├─ ∂ϕodf : ", typeof(cache.∂ϕodf))
    if ~isnothing(cache.∂ϕodf)
        println(io, "      └─── size = ", size(cache.∂ϕodf))
    end
    println(io,   " └─ odf   : ", typeof(cache.odf))
    if ~isnothing(cache.odf)
        println(io, "      └─── size = ", size(cache.odf))
    end
end