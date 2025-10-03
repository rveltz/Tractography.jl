using Tractography
using Test

@testset "Tractography.jl" begin
    include("basic-sampling.jl")
    include("cache.jl")
    include("data.jl")
    include("sh.jl")
end
