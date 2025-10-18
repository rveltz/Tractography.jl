# Structural connectivity estimate

```@contents
Pages = ["connectivity.md"]
Depth = 3
```

The algorithm [`Tractography.Connectivity`](@ref) allows to estimate structural connectivity. In short, 

```julia
alg = Tratography.Connectivity(Tratography.Deterministic())
streamlines, tract_length = TG.sample(model, alg, seeds);
```

returns the first/last point on the streamline. You are then required to map both points with an atlas.