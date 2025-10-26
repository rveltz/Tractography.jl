# Overview of capabilities

A complete streamline generation workflow typically involves the following steps:

1. **Estimating fibre orientation distributions (FODs)**  
2. **Generating seeds** from the estimated FODs  
3. **Running tractography algorithms** to produce streamlines  

The `Tractography.jl` package focuses on steps **2** and **3**.  
Step **1** — estimating FODs — should be performed using a dedicated library, such as [Fibers.jl](https://github.com/lincbrain/Fibers.jl).


---

## FOD specification

The FOD can be supplied to the Tractography Markov Chain (TMC) by the keyword argument `foddata`. See [`Tractography.TMC`](@ref) for more details. Now, you can pass to `foddata` a `Tractography.FODData` object (see [`Tractography.FODData`](@ref) for more information), which can be created:
- from a NIfTI file using the simplified `Tractography.FODData` constructor.
- from an `AbstractArray` object using directly the `Tractography.FODData` constructor.

---

## Seed Generation

Seeds can be generated directly from the FOD or from a binary mask using the following methods:

- [`Tractography.from_fod`](@ref)  
- [`Tractography.from_mask`](@ref)

---

## Tracking Algorithms

- Four tracking algorithms are currently implemented (see [Algorithms](@ref algos)).  
- All algorithms can run on either **GPU** or **CPU (threaded)**, and are **vendor-agnostic**.  
- The interface supports **precomputed caches**, making the algorithms **allocation-free** during execution. This allows for efficient, long-duration runs.  
- Computations support **arbitrary floating-point precision** (`Float32`, `Float64`, `BigFloat`, ...).

---

## Plotting

The package provides convenient plotting utilities to:

- Visualize streamlines  
- Display FODs using glyph-based representations  

---

## Exporting Streamlines

Generated streamlines can be exported in the **`.tck`** format.  
Alternatively, users can export data to any preferred format using standard Julia packages (e.g. `JLD2.jl`).