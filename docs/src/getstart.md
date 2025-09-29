# 🚀 Get started with with Tractography.jl
```@contents
Pages = ["getstart.md"]
Depth = 3
```


This tutorial will introduce you to the functionalities for computing streamlines.

# Basic use

In this example, we will sample `Nmc` streamlines from a Tractography Markov Chain (TMC).

```@example GS
using Tractography
const TG = Tractography

model = TMC(Δt = 0.125f0,
            odfdata = ODFData((@__DIR__) * "/../../examples/fod-FC.nii.gz"),
            )
Nmc = 10
seeds = rand(Float32, 6, Nmc)
alg = Probabilistic()
streamlines, tract_length = sample(model, alg, seeds);
size(streamlines)
```

## Step 1: Define a TMC

We define a Tractography Markov Chain (TMC) model as follows:

```@example GS
model = TMC(Δt = 0.125f0,
            odfdata = ODFData((@__DIR__) * "/../../examples/fod-FC.nii.gz"),
            )
```

## Step 2: Define the seeds

```@example GS
Nmc = 10 # Monte Carlo sample
seeds = rand(Float32, 6, Nmc)
```

## Step 3: Chose a sample algorithm

```@example GS
alg = Probabilistic()
```

## Step 4: Sample the streamlines
```@example GS
streamlines, tract_length = sample(model, alg, seeds);
```

## Optional: export the result for `mrview`

```julia
TG.save_streamlines("tractogram-julia.tck", streamlines)
```

# Optimal use

It is often better to cache some data when computing batches of streamlines. This would be done as follows

```@example GS
model = TMC(Δt = 0.125f0,
            odfdata = ODFData((@__DIR__) * "/../../examples/fod-FC.nii.gz"),
            )
Nmc = 10
seeds = rand(Float32, 6, Nmc)
streamlines = zeros(Float32, 6, 20, Nmc)
tract_length = zeros(UInt32, Nmc)
alg = Probabilistic()
cache = TG.init(model, alg)
# this can be called repeatedly after updating seeds for example
TG.sample!(streamlines, tract_length, model, cache, alg, seeds);
size(streamlines)
```