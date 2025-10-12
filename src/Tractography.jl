module Tractography
    using DocStringExtensions
    using Accessors
    using NIfTI
    using Parameters, StatsBase, LinearAlgebra
    import StaticArrays as SA
    using LoopVectorization
    import FastTransforms

    # sampling method of FODF
    export FibonacciSH, ComputeAllODF
    export ODFData, TMC, Cone, sample_odf, sample, init
    export Probabilistic, Rejection, Deterministic, Diffusion, Connectivity
    export save_streamlines

    # plotting
    export plot_streamlines!, plot_odf, plot_odf!

    include("plot.jl")
    include("tmc.jl")
    include("utils.jl")
    include("data.jl")
    include("sph.jl")
    include("tmccache.jl")
    include("sample.jl")
    include("gpu.jl")
    include("diffusion.jl")
end
