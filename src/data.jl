"""
$(TYPEDEF)

Structure to hold the affine transform from the real world to voxel coordinates.

# Fields
$(TYPEDFIELDS)

# Methods
- see `transform(tf::Transform, x) `
"""
struct Transform{𝒯s, 𝒯t}
    "Forward transform."
    S::𝒯s
    "Inverse transform."
    Sinv::𝒯s
    "Translation."
    T::𝒯t
end
@inline transform(tf::Transform, x) = tf.S * x # for plotting
@inline transform(tf::Transform, x::CartesianIndex{3}) = transform(tf, SA.SVector(Tuple(x)..., 1)) # for plotting
@inline transform(tf::Transform, x::SA.SVector{3}) = transform(tf, SA.SVector(x..., 1)) # for plotting
@inline transform_inv(tf::Transform, x) = tf.Sinv * x

"""
$(TYPEDEF)

Structure to hold ODF data.

# Fields
$(TYPEDFIELDS)

# Methods
- `get_lmax(::ODFData)`
- `size(::ODFData)` return the size of the data
- `get_range(::ODFData)` 
"""
struct ODFData{𝒯, 𝒯d, 𝒯s, 𝒯t}
    "filename from which the (odf) data is read."
    filename::String
    data::𝒯d
    "max l coordinate in of spherical harmonics."
    lmax::Int
    "transform associated with data, see `Transform`."
    transform::Transform{𝒯s, 𝒯t}
    "Are the data normalized? In this case odf[i,j,l,1] ∈ {0,1}."
    normalized::Bool
end

"""
max l coordinate in of spherical harmonics
"""
@inline get_lmax(odfdata::ODFData) = odfdata.lmax
Base.size(odfdata::ODFData) = size(odfdata.data)
_get_array(x::AbstractArray) = x
_get_array(x::NIfTI.NIVolume) = x.raw
_get_array(x::ODFData) = _get_array(x.data)
_my_typeof(x) = typeof(x)
_my_typeof(x::NIfTI.NIVolume) = typeof(x.raw)
is_normalized(odfdata) = odfdata.normalized

function get_range(odfdata::ODFData)
    nx, ny, nz, nt = size(odfdata)
    lx, ly, lz = transform(odfdata, SA.SVector(1, 1, 1))
    rx, ry, rz = transform(odfdata, SA.SVector(nx, ny, nz))
    return sort(LinRange(lx, rx, nx)), 
           sort(LinRange(ly, ry, ny)),
           sort(LinRange(lz, rz, nz))
end

"""
$(SIGNATURES)

Constructor for `ODFData` based on Array data and transform.

## Arguments
- `data::AbstractArray{𝒯, 4}`
"""
function ODFData(file_name, data::AbstractArray{𝒯, 4}, S::𝒯s, T::𝒯t, normalize_it::Bool) where {𝒯, 𝒯s, 𝒯t}
    Sinv = pinv(S)
    lmax = get_lmax_from_odf_length(size(_get_array(data), 4))
    ODFData{_my_typeof(data), typeof(data), 𝒯s, 𝒯t}(file_name, data, lmax, Transform(S, Sinv, T), normalize_it)
end
@inline transform(ni::ODFData, x) = transform(ni.transform, x)
@inline transform_inv(ni::ODFData, x) = transform_inv(ni.transform, x)

"""
$(SIGNATURES)

Constructor for `ODFData` based on NII file.
Read a `.nii.gz` or a `.nii` file passed as a `String`.

The raw spherical harmonics are scaled so that the zero spherical harmonic coefficient is one (or zero).

You can display more information using 

```
show(stdout, ni; full = true)
```

## Output

It returns a `ODFData` struct.
"""
function ODFData(file::String; normalize_it::Bool = true, k...) 
    data = niread(file; k...)
    # we normalize the ODF to have mass one
    if ~all(x-> x >= 0, data.raw[:,:,:,1]) 
        @warn "Some zero SH coefficients are negative!\nPutting them to zero"
    end
    if normalize_it
        _normalize_sph_data!(data)
    end
    @debug size(data) data.header
    A = NIfTI.getaffine(data.header)
    S = SA.@SMatrix [A[i, j] for i = 1:4, j = 1:4]
    T = SA.@SVector [A[i, end] for i = 1:3]
    return ODFData(file, data, S, T, normalize_it)
end

function _normalize_sph_data!(data)
    nx,ny,nz, = size(data.raw)
    Threads.@threads for k=1:nz
        for j=1:ny
            for i=1:nx
                α = data.raw[i,j,k,1]
                if α > 0
                    data.raw[i,j,k,:] ./= α
                end
            end
        end
    end
end
###########################################################################
"""
$(SIGNATURES)
"""
function Base.show(io::IO, ni::ODFData{T, Tp}; full::Bool = false, prefix = "") where {T, Tp} 
    printstyled(prefix, Tp, "\n", bold = true, color = :cyan)
    println(prefix * " ├─ File name   = ", ni.filename)
    println(prefix * " ├─ lmax (SH)   = ", ni.lmax)
    println(prefix * " ├─ Dimensions  = ", size(ni.data))
    println(prefix * " ├─ normalized  = ", ni.normalized)
    if ni.data isa NIfTI.NIVolume
        println(prefix * " ├─ Voxel size  = ", ni.data.header.pixdim[1:4])
        println(prefix * " ├─ Orientation = ", NIfTI.orientation(ni.data))
    end
    println(prefix * " └─ Transform (s_row) = ⋯")
    if full
        ni.transform.S |> display
    end
end