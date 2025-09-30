# GPU example

Here is an example which works on GPU. The implementation is agnostic to the GPU vendor but we focus on NVIDIA for this example.

To speed up computations, we restrict to `Float32`, this is handled easily by the package.

!!! danger
    Because allocations matter a lot on GPU, we pre-allocate all the memory required for computations and fill this inplace.

```@example GPU
using Tractography
const TG = Tractography

# define the model for TMC
model = TMC(Δt = 0.125f0,
            odfdata = ODFData((@__DIR__) * "/../../examples/fod-FC.nii.gz"),
            cone = Cone(45f0),
            proba_min = 0.005f0,
            )
```

# Define the seeds

```julia
using CUDA
# number of streamlines
Nmc = 1024*400
# maximum number of steps for each streamline
Nt = 2000
# define the seeds
seeds = cu(zeros(6, Nmc));
seeds[1:3, :] .= [-13.75, 26.5, 8] .+ 0.1  .* randn(3, Nmc) .|> Float32 |> CuArray;
seeds[4, :] .= 1
tract_length = CuArray((zeros(UInt32, Nmc))
```

we next define a buffer to hold the streamlines

```julia
streamlines_gpu = cu(zeros(Float32, 3, Nt, Nmc), unified = true)
```

# Define the computation cache

Because we can compute the streamlines in batches for the same TMC, it is best to cache some data for this

```julia
# we precompute the cache which is heavy otherwise each call to sample
# will recompute it
cache_g = TG.init(model, Probabilistic(); 
                  𝒯ₐ = CuArray,
                  n_sphere = 400);
```

# Compute the streamlines

The following takes 0.5s on a A100.

```julia
# this setup works for a GPU with 40GiB
# it yields 1e6/sec streamlines for Probabilistic
# and 2.2e6/sec streamlines for Deterministic
CUDA.@time TG.sample!(
            streamlines_gpu,
            tract_length,
            model,
            cache_g,
            Probabilistic(),
            # Deterministic(),
            seeds;
            # maxodf_start = true,
            # reverse_direction = true,
            gputhreads = 1024,
            );
```

!!! tip "batches"
    This can be called many times, for example after updating the seeds.


So far the streamlines `streamlines_gpu` are on the GPU. We can retrieve them at zero cost on the CPU

```julia
streamlines = @time unsafe_wrap(Array, streamlines_gpu);
```