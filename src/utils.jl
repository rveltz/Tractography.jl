####################################################################################################
macro time_debug(msg, ex)
    quote
        local ret = @timed $(esc(ex))
        local _msg = $(esc(msg))
        local _msg_str = _msg === nothing ? _msg : string(_msg)
        time = strip(sprint(Base.time_print, ret.time*1e9, ret.gcstats.allocd, ret.gcstats.total_time, Base.gc_alloc_count(ret.gcstats)))
        @debug _msg_str * " " * time
        ret.value
    end
end
####################################################################################################
"""
$(SIGNATURES)

Transform spherical to cartesian coordinates.

Recall that t∈[0, π] and p∈[0, 2π]
"""
@inline function spherical_to_euclidean(θ, ϕ)
    st, ct = sincos(θ)
    sp, cp = sincos(ϕ)
    x = st * cp
    y = st * sp
    z = ct
    return x, y, z 
end

"""
$(SIGNATURES)

Transform cartesian to spherical coordinates. Assume that the vector has norm one.
"""
@inline function euclidean_to_spherical(x, y, z)
    t = acos(z)
    p = atan(y, x)
    return t, p
end
####################################################################################################
function fibonacci_sampling(N, 𝒯::DataType = Float64)
    out = Vector{Tuple{𝒯, 𝒯}}(undef, N+1)
    ϕ = (1 + √5)/2
    I = 0:N .+ 1/2
    r = 2π/ϕ
    n = 1
    for i in I
        out[n] = (𝒯(acos(1 - 2i/N)), 𝒯(mod(r * i, 2π)))
        n += 1
    end
    out
end
####################################################################################################
@inline softplus(x, k = 1) = ifelse(k * x < 30, log1p(exp(k * x)) / k, x) # avoid Inf for x large
@inline ∂softplus(x, k = 1) = 1 / (1 + exp(-k * x))

"""
Return length of each streamline.
"""
@views function _get_length(streamlines::AbstractArray{𝒯}) where 𝒯
    @assert size(streamlines, 1) >= 3
    _,nt,Nmc = size(streamlines)
    result = zeros(Int, Nmc)
    Threads.@threads for n = 1:Nmc
        len_st = 1
        @inbounds for t = 1:nt-1
            nm = norm(SVector(streamlines[1:3,t,n]...) - SVector(streamlines[1:3,t+1,n]...))
            if nm == 0 || isnan(nm)
                break
            else
                len_st += 1
            end
        end
        result[n] = len_st
    end
    result
end