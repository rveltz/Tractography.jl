abstract type AbstractCache end
# evaluation of spherical harmonics
abstract type AbstractSPHEvaluation end

"""
$(TYPEDEF)

The evaluation of spherical harmonics is done on the fly. Requires little memory.

See also `PreComputeAllODF`
"""
struct DirectSH <: AbstractSPHEvaluation end

"""
$(TYPEDEF)

Set up for plotting ODF.
"""
struct PlottingSH <: AbstractSPHEvaluation end

"""
$(TYPEDEF)

Spherical harmonics evaluation based on Fibonacci sampling. All ODF are pre-computed once and saved in a cache. Their positivity is enforced with a `max(0,⋅)` or a mollifier. 

!!! danger 
    Requires a relatively large memory!

## Details
If you have `na` angles for sampling the unit sphere and the data is of size `(nx, ny, nz, nsph)`, it yields a matrix of dimensions `(nx, ny, nz, na)`.
"""
struct PreComputeAllODF <: AbstractSPHEvaluation end
####################################################################################################
# streamlines tracking algorithms
abstract type AbstractSampler end
# sampler that are based on a grid. Basically everything except ::Rejection
abstract type AbstractNotPureRejectionSampler <: AbstractSampler end
# Deterministic samplers
abstract type DeterministicSampler <: AbstractNotPureRejectionSampler end
abstract type AbstractSDESampler <: AbstractSampler end

"""
$(TYPEDEF)

Tractography sampling performed with the cumulative sum distribution. Can be used with `FibonacciSH` and `PreComputeAllODF`.

# Constructor

`Probabilistic()`
"""
struct Probabilistic <: AbstractNotPureRejectionSampler end

"""
$(TYPEDEF)

Tractography sampling performed with the argmax function. Can be used with `FibonacciSH` and `PreComputeAllODF`.
"""
struct Deterministic <: DeterministicSampler end

"""
$(TYPEDEF)

Tractography based sampling of structural connectivity. 
Do not compute the full streamline but only return the first/last points and the streamlines lengths. This allows to compute many more streamlines on GPU where memory is limited.

## Constructor example
 - `Connectivity(Probabilistic())`
"""
struct Connectivity{Talg} <: AbstractSampler
    alg::Talg
end

_get_alg(alg) = alg
_get_alg(alg::Connectivity) = alg.alg

"""
$(TYPEDEF)

Tractography sampling performed with diffusive model. Basically, the streamlines (Xₜ)ₜ are solution of the SDE

dXₜ = γ * drift(Xₜ) dt + γ_noise * √γ * dnoiseₜ

where

drift(Xₜ) = (Xₜ², ∇log f(Xₜ))

# Arguments (with default values):
$(TYPEDFIELDS)

# Constructor

Example for `Float32`: `Diffusion(γ = 1f0)`.
If you want `Float64`, you have to pass the two scalars

    ```Diffusion(γ = 1.0, γ_noise = 1.0)```
"""
@with_kw_noshow struct Diffusion{Ta, T, Tk, Tmol, Tdmol} <: AbstractSDESampler
    "SciML algorithm used to simulate the tractography diffusion process."
    alg_sde::Ta = nothing
    "γ parameter of the diffusion process."
    γ::T = 1f0
    "parameter of the diffusion process to scale the variance."
    γ_noise::T = 1f0
    "keyword arguments passed to alg_sde."
    kw_sde::Tk = nothing
    "mollifier."
    mollifier::Tmol = Base.Fix2(softplus, 10)
    "differential of mollifier."
    d_mollifier::Tdmol = Base.Fix2(∂softplus, 10)
    "Fixed time step?"
    adaptive::Bool = false
end

get_γ(alg::AbstractSDESampler) = alg.γ
get_γ(alg::Connectivity) = get_γ(_get_alg(alg))
get_γ_noise(alg::AbstractSDESampler) = alg.γ_noise
get_γ_noise(alg::Connectivity) = get_γ_noise(_get_alg(alg))
is_adaptive(alg::AbstractSDESampler) = alg.adaptive
is_adaptive(alg::Connectivity) = is_adaptive(_get_alg(alg))

"""
$(TYPEDSIGNATURES)

Define a transport algorithm. Options are the same as for `Diffusion`.
"""
function Transport(;kwargs...)
    alg = Diffusion(;kwargs...)
    @reset alg.γ_noise = zero(alg.γ_noise)
end

function Base.show(io::IO, alg::Diffusion{Ta, T}) where {Ta, T}
    printstyled(io, "Diffusion [$T]" ; bold = true, color = :cyan)
    printstyled(io, " sampling algorithm\n", color = :cyan)
    println(io, "├─ adaptive = ", is_adaptive(alg))
    println(io, "├─ alg      = ", alg.alg_sde)
    println(io, "├─ γ        = ", get_γ(alg))
    println(io, "└─ γ_noise  = ", get_γ_noise(alg))
    if ~isnothing(alg.kw_sde)
        println(io, "kw = ", alg.kw_sde)
    end
end
####################################################################################################
"""
$(TYPEDEF)

Structure to encode a cone to limit sampling the direction. 
This ensures that the angle in degrees between to consecutive streamline directions is less than `angle`.

The implemented condition is for `cn = Cone(angle)`.

```
(cn::Cone)(d1, d2) = dot(d1, d2) > cosd(cn.alpha)
```

## Fields

$(TYPEDFIELDS)

## Constructor

`Cone(angle)`
"""
struct Cone{𝒯 <: Real}
    "half section angle in degrees"
    alpha::𝒯
end
(cn::Cone)(d1, d2) = dot(d1, d2) > cosd(cn.alpha)

"""
$(TYPEDEF)

Tractography Markov Chain (TMC).

# Fields (with default values):
$(TYPEDFIELDS)

# Methods
- `_apply_mask!(model, mask)` apply a mask to the raw SH tensor. See its doc string.
- `_getdata(model)` return the fodf data associated with the TMC.
- `size(model)` return `nx, ny, nz, nt`.
- `eltype(model)` return the scalar type of the data (default Float64).
- `get_lmax(model)` return the max `l` coordinate in of spherical harmonics.

# Constructors (use the fields!)
- `TMC()`
- `TMC(Δt = 0.1f0)` for a Float32 TMC
- `TMC(Δt = 0.1, proba_min = 0.)` for a Float64 TMC. You need to specify both fields `Δt` and `proba_min`
- `TMC(odfdata = rand(10,10,10,45))` for custom ODF
"""
@with_kw_noshow struct TMC{𝒯, 𝒯alg <: AbstractSPHEvaluation, 𝒯d, 𝒯C, 𝒯mol}
    "Step size of the TMC."
    Δt::𝒯 = 0.1f0
    "Spherical harmonics evaluation algorithm. Can be `FibonacciSH(), PreComputeAllODF()`."
    evaluation_algo::𝒯alg = PreComputeAllODF()
    "ODF data from nifti file. Must be the list of ODF in the base of spherical harmonics. Hence, it should be an (abstract) 4d array."
    odfdata::𝒯d = nothing
    "Cone function to restrict angle diffusion. You can use a `Cone` or a custom function `(d1, d2) -> return_a_boolean`."
    cone::𝒯C = Cone(90f0)
    "Probability below which we stop tracking."
    proba_min::𝒯 = 0.0f0
    "Mollifier, used to make the fodf non negative. During odf evaluation, we effectively use `mollifier(fodf[angle,i,j,k])`."
    mollifier::𝒯mol = max_mollifier
end
@inline getdata(model::TMC) = model.odfdata
Base.size(model::TMC) = size(getdata(model))
Base.eltype(model::TMC{𝒯}) where 𝒯 = 𝒯
@inline get_lmax(model::TMC) = get_lmax(getdata(model))
"""
$(SIGNATURES)

`max(x,0)` as mollifier to prevent negative ODF.
"""
max_mollifier(x) = max(0, x)
get_range(model::TMC) = get_range(getdata(model))
get_array(model::TMC) = _get_array(getdata(model))

function Base.show(io::IO, model::TMC)
    printstyled(io, "TMC with elype ", eltype(model), bold = true, color = :cyan)
    println(io, "\n ├─ Δt = ", model.Δt)
    println(io, " ├─ minimal probability = ", model.proba_min)
    if model.cone isa Cone
        println(io, " ├─ cone                = ", model.cone)
    end
    println(io, " ├─ mollifier           = ", model.mollifier)
    println(io, " ├─ evaluation of SH    = ", model.evaluation_algo)
    if model.odfdata isa ODFData
        println(io, " └─ data : (lmax = $(get_lmax(model)))")
        show(io, model.odfdata; prefix = "      ")
    end
    if model.odfdata isa AbstractArray
        println(io, " └─ data                = ", typeof(model.odfdata))
    end
end

"""
$(SIGNATURES)

Multiply the mask which is akin to a matrix of `Bool` with same size as the data stored in `model`. Basically, the mask `mask[ix, iy, iz]` ensures whether the voxel `(ix, iy, iz)` is discarded or not.

# Arguments

- `model::TMC`.
- `mask` can be a `AbstractArray{3, Bool}` or a `NIVolume`.
"""
function _apply_mask!(model, mask)
    if ~isnothing(mask)
        nx, ny, nz, nsh = size(model)
        data = _get_array(getdata(model))
        for k = 1:nsh
            @tturbo data[:, :, :, k] .*= mask
        end
    end
end
####################################################################################################
function save_streamlines end