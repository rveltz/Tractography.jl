abstract type AbstractCache end
# evaluation of Spherical harmonics
abstract type AbstractEvaluation end

"""
$(TYPEDEF)

Spherical harmonics evaluation based on Fibonacci sampling. 
The ODF are computed on the fly using matrix-vector. Requires little memory.

See also `ComputeAllODF`
"""
struct FibonacciSH <: AbstractEvaluation end

"""
$(TYPEDEF)

Spherical harmonics evaluation based on Fibonacci sampling. All ODF are pre-computed once. Their positivity is enforced with a `max(0,⋅)` or a mollifier. 

!!! danger 
    Requires a relatively large memory!

## Details
If you have `na` angles for sampling the unit sphere and the data is of size `(nx,ny,nz,nsph)`, it yields a matrix of dimensions `(nx,ny,nz,na)`.

See also `FibonacciSH`
"""
struct ComputeAllODF <: AbstractEvaluation end
####################################################################################################
# streamlines tracking algorithms
abstract type AbstractSampler end
# sampler that are based on a grid. Basically everything except ::Rejection
abstract type AbstractNotPureRejectionSampler <: AbstractSampler end
# Deterministic samplers
abstract type DeterministicSampler <: AbstractNotPureRejectionSampler end

"""
$(TYPEDEF)

Tractography sampling performed with the cumulative sum distribution. Can be used with `FibonacciSH` and `ComputeAllODF`.

# Constructor

`CSD()`
"""
struct CSD <: AbstractNotPureRejectionSampler end

"""
$(TYPEDEF)

Tractography sampling performed with the argmax function. Can be used with `FibonacciSH` and `ComputeAllODF`.
"""
struct Deterministic <: DeterministicSampler end

"""
$(TYPEDEF)

Tractography based sampling of Connectivity. 
Do not compute the full streamline, only return first/last points and streamlines lengths.

## Constructor example
 - `Connectivity(CSD())`
"""
struct Connectivity{Talg} <: AbstractSampler
    alg::Talg
end

_get_alg(alg) = alg
_get_alg(alg::Connectivity) = alg.alg
####################################################################################################
"""
$(TYPEDEF)

Structure to encode a cone to limit sampling the direction. 
This ensures that the angle in degrees between to consecutive streamline directions is less than `angle`.

The implemented condition is for `cn = Cone(angle)`

```
(cn::Cone)(d1, d2) = dot(d1, d2) > cos(pi/180 * cn.alpha)
```

## Fields

$(TYPEDFIELDS)

## Constructor

`Cone(angle)`
"""
struct Cone{T <: Real}
    "half section angle in degrees"
    alpha::T
end
(cn::Cone)(d1, d2) = dot(d1, d2) > cos(pi/180 * cn.alpha)

"""
$(TYPEDEF)

Tractography Markov Chain (TMC).

# Fields (with default values):
$(TYPEDFIELDS)

# Methods
- `_apply_mask!(model, mask)` apply a mask to the raw SH tensor. See its doc string.
- `_getdata(model)` return the fodf data associated with the TMC.
- `size(model)` return `nx,ny,nz,nt`.
- `eltype(model)` return the scalar type of the data (default Float64).
- `get_lmax(model)` return the max `l` coordinate in of spherical harmonics.

# Constructors (use the fields!)
- `TMC()`
- `TMC(Δt = 0.1f0)` for a Float32 TMC
- `TMC(Δt = 0.1, proba_min = 0.)` for a Float64 TMC. You need to specify both fields `Δt` and `proba_min`
- `TMC(odfdata = rand(10,10,10,45))` for custom ODF
"""
@with_kw_noshow struct TMC{T, Talg <: AbstractEvaluation, Td, TC, Tmol}
    "Step size of the TMC."
    Δt::T = 0.1f0
    "Spherical harmonics evaluation algorithm. Can be `FibonacciSH(), ComputeAllODF()`."
    evaluation_algo::Talg = ComputeAllODF()
    "ODF data from nifti file. Must be the list of ODF in the base of spherical harmonics. Hence, it should be an (abstract) 4d array."
    odfdata::Td = nothing
    "Cone function to restrict angle diffusion. You can use a `Cone` or a custom function `(d1,d2) -> return_a_boolean`."
    C::TC = Cone(90)
    "Probability below which we stop tracking."
    proba_min::T = 0.0f0
    "Mollifier, used to make the fodf non negative. During odf evaluation, we effectively use `mollifier(fodf[i,j,k,angle])`."
    mollifier::Tmol = default_mollifier
end
@inline getdata(model::TMC) = model.odfdata
Base.size(model::TMC) = size(getdata(model))
Base.eltype(model::TMC{T}) where T = T
@inline get_lmax(model::TMC) = get_lmax(getdata(model))
default_mollifier(x) = max(0, x)
get_range(model::TMC) = get_range(getdata(model))
get_array(model::TMC) = _get_array(getdata(model))

function Base.show(io::IO, model::TMC)
    printstyled(io, "TMC with elype ", eltype(model), bold = true, color = :cyan)
    println(io, "\n ├─ Δt = ", model.Δt)
    println(io, " ├─ minimal probability = ", model.proba_min)
    if model.C isa Cone
        println(io, " ├─ cone                = ", model.C)
    end
    if model isa TMC
        println(io, " ├─ mollifier           = ", model.mollifier)
    end
    if model isa TMC
        println(io, " ├─ evaluation of SH    = ", model.evaluation_algo)
        if model.odfdata isa ODFData
            println(io, " └─ data : ⋯")
            show(io, model.odfdata; prefix = "      ")
        end
        if model.odfdata isa AbstractArray
            println(io, " └─ data                = ", typeof(model.odfdata))
        end
    end
end

"""
$(SIGNATURES)

Multiply the mask which is akin to a matrix of `Bool` with same size as the data stored in `model`. Basically, the mask `mask[ix,iy,iz]` ensures whether the voxel `(ix,iy,iz)` is discarded or not.

# Arguments

- `model::TMC`.
- `mask` can be a `AbstractArray{3, Bool}` or a `NIVolume`.
"""
function _apply_mask!(model, mask)
    if ~isnothing(mask)
        nx, ny, nz, nsh = size(model)
        data = _get_array(getdata(model))
        for k = 1:nsh
            @tturbo data[:,:,:,k] .*= mask
        end
    end
end
####################################################################################################
"""
$(SIGNATURES)

Save tractogram in `filename` in a `tck` format. You should probably use a filename like `streamlines-julia.tck`.

# Arguments
- `filename::String` file name.
- `tractogram` can either be an abstract array of dimension `3 x n_length_fiber x n_fibres` or a `Vector{ <: TractoResult}` after a call to `sample`.

# Notes
The current implementation is based on NiBabel.
"""
function save_streamlines(filename::String, tractogram::AbstractArray{T, 3}, tracto_length) where {T}
    if size(tractogram, 1) != 3 
        error("The streamlines must have size 3 x ns x n_nmc")
    end
    nib = pyimport("nibabel")
    np = pyimport("numpy")
    streamlines = tracto = [tractogram[:, 1:tracto_length[k], k]' for k in axes(tractogram,3) if tracto_length[k]>2]
    @info "Saving in file: " filename
    nib.streamlines.TckFile(nib.streamlines.Tractogram(streamlines, affine_to_rasmm = np.eye(4))).save(filename)
end