# Overview of capabilities

A general streamline generation pipeline consists of 

1. estimate fibre orientation distributions (FOD),
2. generate seeds from the FOD,
3. run the tractography algorithm to generate streamlines.

The current package `Tractography.jl` addresses points 2 and 3, step 1 has to be performed by another library, for example [Fibers.jl](https://github.com/lincbrain/Fibers.jl).


## Seeds generation

The seeds can be generated using the methods [`Tractography.from_odf`](@ref) and [`Tractography.from_mask`](@ref). 

## Tracking algorithms

- 4 tracking algorithms are provided, see [pages algos](@ref algos).
- the algorithms run on GPU or on CPU (threaded), vendor agnostic.
- interface is provided such that once a cache has been precomputed, the algorithms are allocation free. They can thus be ran as long as necessary.
- the computation can be run in arbitrary precision `Float32, Float64, BigFloat, ...`. 

## Plotting

We provide plotting capabilities to 
- view streamlines,
- plot the FODs using glyphs. 

## Exporting streamlines

We provide a method to export the results in `tck` format. Of course, the user can use Julia package (*e.g.*  `JLD2.jl`) to export to their preferred format.