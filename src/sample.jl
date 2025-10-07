function _init_fibonacci_sh(model::TMC{𝒯}, n_sphere) where {𝒯}
    lmax = get_lmax(model)
    angles = fibonacci_sampling(n_sphere, 𝒯)
    directions = [spherical_to_euclidean(d[1], d[2]) for d in angles]
    n_angles = length(angles)
    Yₗₘ = get_vector_of_sh(angles, lmax)
    cone = isnothing(model.cone) ? nothing : [𝒯(model.cone(d1, d2)) for d1 in directions, d2 in directions]
    TMCCache(; n_sphere = n_angles, Yₗₘ, dΩ = 𝒯(4pi / n_angles), angles, lmax, cone, directions)
end

# for plotting
function _init(model::TMC{𝒯, PlottingSH},
                alg::AbstractNotPureRejectionSampler; 
                n_sphere = 400) where 𝒯
    _init_fibonacci_sh(model, n_sphere)
end

function _init(model::TMC{𝒯, PreComputeAllODF},
                alg::AbstractNotPureRejectionSampler; 
                n_sphere = 400) where 𝒯
    cache = _init_fibonacci_sh(model, n_sphere)
    _build_cache_from_Y_matrix(model, cache, cache.Yₗₘ)
end

function _build_cache_from_Y_matrix(model::TMC{𝒯}, cache, Ysv) where 𝒯
    na = size(Ysv, 1)
    # compute all ODF
    nx, ny, nz, nt = size(model)
    # we optimize the array layout because loops occurs on the angle variable
    # and julia is column major
    ni = 𝒯.(get_array(model))
    ni_v = reshape(ni, nx*ny*nz, nt)
    odf_v  = @time_debug "Mat-Vec:" ni_v * Ysv'; # vector view
    odf = reshape(odf_v, nx, ny, nz, na);
    @time_debug "Mollifier:" @tturbo @. odf = model.mollifier(odf)
    @reset cache.odf = permutedims(odf, (4, 1, 2, 3))
    return cache
end

"""
$(TYPEDSIGNATURES)

Sample the TMC `model`.

## Arguments
- `model::TMC`
- `alg` sampling algorithm, `Deterministic, Probabilistic, Diffusion, etc`.
- `seeds` matrix of size `6 x Nmc` where `Nmc` is the number of Monte-Carlo simulations to be performed.
- `mask = nothing` matrix of boolean where to stop computation. See also `_apply_mask`.

## Optional arguments
- `nt::Int` maximal length of each streamline.
- `n_sphere::Int = 400` number of points to discretize the sphere on which we evaluate the spherical harmonics.
- `maxodf_start::Bool` for each locations, use direction provided by the argmax of the ODF.
- `reverse_direction::Bool` reverse initial direction.
- `nthreads::Int = 8` number of threads on CPU.
- `gputhreads::Int = 512` number of threads on GPU.

## Output
- `streamlines` with shape `3 x nt x Nmc`
"""
function sample(model::TMC{𝒯},
                alg::AbstractSampler,
                seeds::AbstractMatrix{𝒯},
                mask::Union{Nothing, AbstractArray{Bool}} = nothing;
                nt::Int = 1000,
                n_sphere::Int = 400,
                maxodf_start::Bool = false,
                reverse_direction::Bool = false,
                nthreads::Int = 8,
                gputhreads::Int = 512
                ) where {𝒯}
    @assert size(seeds, 1) == 6
    cache = init(model, alg; n_sphere)
    if ~isnothing(mask)
        _apply_mask!(model, mask)
    end
    streamlines = similar(seeds, 3, (alg isa Connectivity ? 2 : nt), size(seeds, 2))
    streamlines_length = zeros(UInt32, size(seeds, 2))
    sample!(streamlines,
            streamlines_length,
            model,
            cache,
            alg,
            seeds;
            maxodf_start,
            reverse_direction,
            nthreads,
            gputhreads,
            nₜ = nt)
    return streamlines, streamlines_length
end