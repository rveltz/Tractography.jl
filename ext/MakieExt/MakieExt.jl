module MakieExt
    using Makie, Tractography, DocStringExtensions, LinearAlgebra, StaticArrays, Accessors
    import Tractography: add_frame!,
                           plot_streamlines!, 
                           plot_odf!,
                           plot_odf,
                           plot_slice,
                           plot_slice!,
                           transform,
                           init, _init,
                           TMC,
                           getdata,
                           _get_sphere_fibonacci,
                           analyse_streamline,
                           get_voxel_gpu, euclidean_to_spherical,
                           spherical_to_euclidean,
                           @time_debug
    include("plot.jl")
end
