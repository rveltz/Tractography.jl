"""
$(SIGNATURES)

Exponential map on the sphere.

We assume that t>0

See https://github.com/JuliaManifolds/ManifoldsBase.jl/blob/5c4a61ed3e5e44755a22f7872cb296a621905f87/test/ManifoldsBaseTestUtils.jl#L63
"""
function Exp𝕊²(p, X, t)
    n = norm(X)
    if iszero(n)
        return p
    end
    s, c = sincos(t * n)
    return c .* p .+ X .* (s / n)
end
####################################################################################################
function init(model::TMC{𝒯, DirectSH}, 
                alg::Union{Talg, Connectivity{ Talg}}; 
                𝒯ₐ = Array{𝒯},
                n_sphere = 0) where {𝒯, Talg <: AbstractSDESampler}
    ni =  𝒯.(get_array(model))
    cache_cpu = TMCCache(; n_sphere, angles = 0, lmax = get_lmax(model), dΩ = zero(𝒯))
    fod = permutedims(ni, (4, 1, 2, 3))
    if is_normalized(model.foddata)
        fod ./= 𝒯(sqrt(4pi))
    end

    ThreadedCache(
        𝒯ₐ(fod),
        𝒯ₐ(zeros(𝒯,0,0,0,0)),
        𝒯ₐ(zeros(𝒯,0,0,0,0)),
        nothing,
        𝒯ₐ(zeros(𝒯,0,0)),
        nothing,
        𝒯ₐ(zeros(𝒯,0,0,0)),
        nothing,
    )
end

function _init(model::TMC{𝒯, PreComputeAllFOD}, 
                alg::AbstractSDESampler; 
                n_sphere = 400) where 𝒯
    # we want to differentiate wrt (θ,ϕ) the expression mollifier(fodf(θ,ϕ))
    # the expression is ∂mollifier(fodf(θ,ϕ)) * ∂fodf(θ,ϕ)
    mollifier = alg.mollifier
    d_mollifier = alg.d_mollifier

    cache = _init_fibonacci_sh(model, n_sphere)
    angles = cache.angles
    na = n_sphere + 1
    lmax = get_lmax(model)
    ∂θYₗₘ = get_vector_of_sh(angles, lmax, 1)
    ∂ϕYₗₘ = get_vector_of_sh(angles, lmax, 2)

    # compute all FOD
    nx, ny, nz, nt = size(model)
    ni =  get_array(model)
    ni_v = 𝒯.(reshape(ni, nx*ny*nz, nt)) # vector version
    Y_v = cache.Yₗₘ
    odf_v  = @time_debug "Mat-Vec:" ni_v * Y_v';
    odf = reshape(odf_v, nx, ny, nz, na);

    # compute all ∂θODF
    ∂θY = 𝒯.(∂θYₗₘ)
    odf_vt = @time_debug "all ∂θodf:" ni_v * ∂θY';
    d_mollifier_odf_v = LV.@tturbo @. d_mollifier(odf_v)
    @time_debug "Apply mollifier" LV.@tturbo @. odf_vt = d_mollifier_odf_v * odf_vt
    ∂θodf = reshape(odf_vt, nx, ny, nz, na);

    # compute all ∂ϕODF
    ∂ϕY = 𝒯.(∂ϕYₗₘ)
    odf_vp = @time_debug "all ∂ϕodf:" ni_v * ∂ϕY';
    @time_debug "Apply mollifier" LV.@tturbo @. odf_vp = d_mollifier_odf_v * odf_vp
    ∂ϕodf = reshape(odf_vp, nx, ny, nz, na);

    @time_debug "Apply mollifier" LV.@tturbo @. odf = mollifier(odf)
    @reset cache.odf   = @time_debug"permutedims" permutedims(odf,   (4, 1, 2, 3))
    @reset cache.∂θodf = permutedims(∂θodf, (4, 1, 2, 3))
    @reset cache.∂ϕodf = permutedims(∂ϕodf, (4, 1, 2, 3))
    @reset cache.∂θYₗₘ = ∂θYₗₘ
    @reset cache.∂ϕYₗₘ = ∂ϕYₗₘ
    return cache
end

function init(model::TMC{𝒯},
                alg::Union{Talg, Connectivity{ Talg}}; 
                n_sphere = 400,
                𝒯ₐ = Array{𝒯},
                ) where {𝒯, Talg <: AbstractSDESampler}
    cache_cpu = _init(model, _get_alg(alg); n_sphere)
    # do not copy the array if the types are the same
    _is_on_cpu = cache_cpu.odf isa 𝒯ₐ
    ∫odf = sum(cache_cpu.odf, dims = 1)[1, :, :, :]
    # here, we have to be careful because the mollifier attributes non zero probabilities
    map!(x -> x > 0 ? x : zero(x), ∫odf, @views model.foddata.data.raw[:,:,:,1])

    ThreadedCache(
            _is_on_cpu ? cache_cpu.odf   : 𝒯ₐ(cache_cpu.odf),
            _is_on_cpu ? cache_cpu.∂θodf : 𝒯ₐ(cache_cpu.∂θodf),
            _is_on_cpu ? cache_cpu.∂ϕodf : 𝒯ₐ(cache_cpu.∂ϕodf),
            _is_on_cpu ? cache_cpu.cone  : 𝒯ₐ(cache_cpu.cone),
            𝒯ₐ(mapreduce(x->[x[1] x[2] x[3]], vcat, cache_cpu.directions)),
            𝒯ₐ(mapreduce(x->[x[1] x[2]],      vcat, cache_cpu.angles)),
            _is_on_cpu ? ∫odf : 𝒯ₐ(∫odf),
            cache_cpu.dΩ
    )
end
####################################################################################################
function sample!(streamlines,
                streamlines_length::AbstractArray{UInt32, 1},
                model::TMC{𝒯}, 
                cache::AbstractCache, 
                alg::Union{AbstractSDESampler, Connectivity{ <: AbstractSDESampler}},
                seeds;
                maxfod_start::Bool = false,
                reverse_direction::Bool = false,
                nthreads = 8,
                gputhreads = 512,
                nₜ = size(streamlines, 2),
                saveat::Int = 1,
                𝒯ₐ = Array) where {𝒯}
    Nmc = size(seeds, 2)
    if size(seeds, 1) != 6 
        error("The initial positions must be passed as an 6 x N array.")
    end
    if (Nmc > size(streamlines, 3))
        error("$Nmc <= ", size(streamlines, 3))
    end
    if ndims(streamlines) < 3
        error("streamlines must be passed as an 3 x nt x N array")
    end
    @debug size(streamlines) nₜ Nmc nthreads gputhreads model.Δt alg

    _, nx, ny, nz = size(cache.odf)
    streamlines_length .= nₜ ÷ saveat
    # launch gpu kernel
    backend = KA.get_backend(seeds)
    nth = backend isa KA.GPU ? gputhreads : nthreads
    kernel! = _sample_kernel_diffusion!(backend, nth)
    @time "kernel-diffusion" kernel!(
                            streamlines,
                            streamlines_length,
                            _get_alg(alg),
                            seeds,
                            cache.odf,
                            cache.∂θodf,
                            cache.∂ϕodf,
                            cache.∫odf,
                            cache.directions,
                            model.foddata.transform,
                            Int32(nₜ),
                            maxfod_start,
                            reverse_direction,
                            model.proba_min,
                            abs(model.Δt),
                            saveat,
                            get_γ(alg),
                            get_γ_noise(alg),
                            cache.dΩ,
                            nx, ny, nz,
                            Val(model.evaluation_algo isa PreComputeAllFOD),
                            Val(~(alg isa Connectivity)),
                            ndrange = Nmc
                            )
    return streamlines
end

KA.@kernel inbounds=false function _sample_kernel_diffusion!(
                            streamlines::AbstractArray{𝒯, 3},
                            streamlines_length::AbstractArray{UInt32, 1},
                            alg::AbstractSDESampler{𝒯},
                            @Const(seeds::AbstractMatrix{𝒯}),
                            @Const( fodf::AbstractArray{𝒯, 4}),
                            @Const(∂θodf::AbstractArray{𝒯, 4}),
                            @Const(∂ϕodf::AbstractArray{𝒯, 4}),
                            @Const(∫odf::AbstractArray{𝒯, 3}),
                            @Const(directions::AbstractMatrix{𝒯}),
                            @Const(tf),
                            @Const(nₜ),
                            @Const(maxfod_start),
                            @Const(reverse_direction),
                            @Const(proba_min::𝒯),
                            @Const(dt::𝒯),
                            @Const(saveat::Int),
                            @Const(γ::𝒯),
                            @Const(γn::𝒯),
                            @Const(dΩ),
                            nx, ny, nz,
                            _precomputed_odf::Val{precomputed_odf},
                            ::Val{save_full_streamline},
                            ) where {𝒯, save_full_streamline, precomputed_odf}
    # index of the streamline being computed
    nₙₘ = @index(Global)

    x₁ = seeds[1, nₙₘ]
    x₂ = seeds[2, nₙₘ]
    x₃ = seeds[3, nₙₘ]
    u₁ = seeds[4, nₙₘ]
    u₂ = seeds[5, nₙₘ]
    u₃ = seeds[6, nₙₘ]

    n_angles::UInt32 = size(directions, 1)

    # current index of angle
    ind_u::Int32 = 1
    # streamline length
    t_length::UInt32 = 1

    if maxfod_start && precomputed_odf
        voxel_index₁, voxel_index₂, voxel_index₃ = get_voxel_index(tf, (x₁, x₂, x₃))
        ind_u = _device_argmax(fodf, voxel_index₁, voxel_index₂, voxel_index₃, n_angles)
        u₁ = directions[ind_u, 1]
        u₂ = directions[ind_u, 2]
        u₃ = directions[ind_u, 3]
    end

    if reverse_direction
        u₁ = -u₁
        u₂ = -u₂
        u₃ = -u₃
    end

    if (reverse_direction || ~maxfod_start) && precomputed_odf
        ind_u = _device_get_angle(directions, u₁, u₂, u₃, n_angles)
    end

    inside_image::Bool = true
    continue_tracking::Bool = true
    voxel_index₁ = voxel_index₂ = voxel_index₃ = Int32(0)

    streamlines[1, 1, nₙₘ] = x₁
    streamlines[2, 1, nₙₘ] = x₂
    streamlines[3, 1, nₙₘ] = x₃

    conditioned_proba = zero(𝒯)
    F = ∫F = Fθ = Fϕ = hx = ∂ = zero(𝒯)
    st = ct = sp = cp = zero(𝒯)
    θᵢ, ϕᵢ = euclidean_to_spherical(u₁, u₂, u₃)
    iₛₐᵥₑ = one(UInt32)

    # Riemannian Langevin algorithm [1]
    # Bharath, Karthik, Karthik Bharath, Alexander Lewis, et al. Sampling and Estimation on Manifolds Using the Langevin Diﬀusion. n.d.
    # X_{n+1}^h =\exp_{X_n^h}(  h/2⋅∇ E(X_n^h) + √h ⋅ g^{-1/2}(X_n^h) ⋅ ξ_{n+1})


    for iₜ = UInt32(2):nₜ
        P = SA.SVector(x₁, x₂, x₃)
        D = SA.SVector(u₁, u₂, u₃)
        # x is in native space
        (voxel_index₁, voxel_index₂, voxel_index₃) = get_voxel_index(tf, (x₁, x₂, x₃))

        inside_image = in_image(voxel_index₁, voxel_index₂, voxel_index₃, nx, ny, nz)
        continue_tracking = inside_image && continue_tracking
        t_length += continue_tracking

        if continue_tracking
            if precomputed_odf
                ind_u = _device_get_angle(directions, u₁, u₂, u₃, n_angles)
                # !! Careful here, we need to have a probability: F / ∫F
                F  =  fodf[ind_u, voxel_index₁, voxel_index₂, voxel_index₃]
                Fθ = ∂θodf[ind_u, voxel_index₁, voxel_index₂, voxel_index₃]
                Fϕ = ∂ϕodf[ind_u, voxel_index₁, voxel_index₂, voxel_index₃]
                ∫F =  ∫odf[voxel_index₁, voxel_index₂, voxel_index₃]
                st, ct = sincos(θᵢ)
                sp, cp = sincos(ϕᵢ)
            else
                F, Fϕ, Fθ = ishtmtx_dot(ϕᵢ, θᵢ, @view fodf[:, voxel_index₁, voxel_index₂, voxel_index₃])
                ∂ = ∂softplus(F, 100f0)
                F =  softplus(F, 100f0)
                Fθ *= ∂
                Fϕ *= ∂
                ∫F = fodf[1, voxel_index₁, voxel_index₂, voxel_index₃]
            end
            continue_tracking = ∫F > proba_min # recall ∫F ∈ [0, 1]
        end

        if continue_tracking
            st, ct = sincos(θᵢ)
            sp, cp = sincos(ϕᵢ)
            # tangent vectors in polar coordinates
            # recall D = (st * cp, st * sp, ct), error ~ 1e-7

            eθ = SA.SVector(ct * cp, ct * sp, -st )
            eϕ = SA.SVector(-sp, cp, 0) # remove the sin(θ) from eϕ because we removed it in Fϕ

            drift = Fθ * eθ + Fϕ * eϕ

            # 19-AAP1507
            if is_adaptive(alg)
                hx = dt * 2 / min(max(1f0, norm(drift)/F)^2, 10)
            else
                hx = dt
            end

            if alg isa Transport
                tangent = (γ * hx / F) * drift
            else
                noise = randn(𝒯) * eθ + randn(𝒯) * eϕ
                tangent = (γ * hx / F) * drift + sqrt(2γ * hx * γn) * noise
            end

            # Geometric-Euler scheme
            u₁, u₂, u₃ = Exp𝕊²(D, tangent, one(𝒯)) # injectivity radius

            # u₁, u₂, u₃ = D + tangent
            # u₁, u₂, u₃ = normalize(SA.SVector(u₁, u₂, u₃))
            θᵢ, ϕᵢ = euclidean_to_spherical(u₁, u₂, u₃)

            x₁ += hx * u₁
            x₂ += hx * u₂
            x₃ += hx * u₃
        else
            streamlines_length[nₙₘ] = t_length ÷ saveat
            if ~save_full_streamline
                streamlines[1, 2, nₙₘ] = x₁
                streamlines[2, 2, nₙₘ] = x₂
                streamlines[3, 2, nₙₘ] = x₃
            end
        end

        if mod(iₜ, saveat) == 0
            if save_full_streamline
                streamlines[1, iₛₐᵥₑ, nₙₘ] = x₁
                streamlines[2, iₛₐᵥₑ, nₙₘ] = x₂
                streamlines[3, iₛₐᵥₑ, nₙₘ] = x₃
                iₛₐᵥₑ += 1
            end
        end
    end
end