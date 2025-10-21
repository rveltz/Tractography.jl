# Plotting 

Plotting is provided by calling recipes to `Makie.jl`. Hence by loading `CairoMakie` or `GLMakie` or `WGLMakie`, you can use custom function to interact with FOD data.
For this part, we use `CairoMakie` but the user is encouraged to use `GLMakie` instead.

Let us start by reading some `nii` data.

```@example PLOTTING
using Tractography
const TG = Tractography

model = TG.TMC(Δt = 0.125f0,
            odfdata = TG.ODFData((@__DIR__) * "/../../examples/fod-FC.nii.gz"),
            cone = TG.Cone(15f0),
            proba_min = 0.005f0,
            )
```

## Plotting the FODFs

We rely on the function [`plot_odf!`](@ref) provides the information regarding the keyword arguments.

```@example PLOTTING
using CairoMakie
f, sc = TG.plot_odf(model; n_sphere = 1500, radius = 1., I = 29:31, J = 30:31, K = 2:2);
cam3d = Makie.cameracontrols(sc) # hide
cam3d.eyeposition[] = Vec3f(96.91587149289712, 85.36944962225634, -85.12319118529459)  # hide
cam3d.lookat[] = Vec3f(95.90360332812566, 85.36944962225634, 2.943478780350152)  # hide
cam3d.upvector[] = Vec3f(0.9999339466437674, 0.0, 0.01149357861675617)  # hide
Makie.zoom!(sc.scene, cam3d, 0.10f0)  # hide
f
```

## Plotting the streamlines

We rely on the function [`plot_streamlines!`](@ref) provides the information regarding the keyword arguments.

```@example PLOTTING
# make some streamlines
streamlines = zeros(Float32, 6, 100, 20)
for n = 1:20
    v0 = rand(3)*15
    for nt = 1:100
        streamlines[1:3, nt, n] .= v0 .+ nt .* [1,1 + 0.1*rand(), 0]
    end
end
# plot the streamlines and the glyph
f = Figure(backgroundcolor = :white)
lscene = LScene(f[1,1])
TG.plot_streamlines!(lscene.scene, streamlines[1:3, 1:1:end, :])
TG.plot_odf!(lscene, model; n_sphere = 1500, radius = 1.3, st = 4);
f
```