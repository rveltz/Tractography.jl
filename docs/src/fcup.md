# FCUP

```@contents
Pages = ["fcup.md"]
Depth = 3
```

This is a more advanced tutorial because we want to show how to apply a mask.


# Define the TMC

```@example FCUP
import Tractography as TG

model = TG.TMC(Δt = 0.125f0,
            odfdata = TG.ODFData((@__DIR__) * "/../../examples/fod-FC.nii.gz"),
            C = TG.Cone(15),
            proba_min = 0.005f0,
            )
```

Just for fun, we plot the FODF of the model.

```@example FCUP
using CairoMakie

f, sc = TG.plot_odf(model; n_sphere = 1500, radius = 0.3, st = 2);
cam3d = Makie.cameracontrols(sc)
cam3d.eyeposition[] = Vec3f(85, 95, -28)
cam3d.lookat[] = Vec3f(84, 95, 59)
rotate_cam!(sc.scene, 0, 0, -pi/2)
f
```

# Define the seeds

We next apply a mask on the boundary of which the streamlines stop.

```@example FCUP
using NIfTI
mask = NIfTI.niread((@__DIR__) * "/../../examples/wm-FC.nii.gz");
TG._apply_mask!(model, mask);
```

We compute `Nmc` streamlines, hence we need `Nmc` seeds

```@example FCUP
using LinearAlgebra

Nmc = 100_000
_ind1 = findall(mask .== 1)
seed_id = rand(1:length(_ind1), Nmc)
seeds = zeros(Float32, 6, Nmc)
for i=1:Nmc
    seeds[:,i] .= vcat(TG.transform(model.odfdata, _ind1[seed_id[i]])[1:3]..., normalize(randn(3)))
end
```

# Compute the streamlines

```@example FCUP
streamlines, tract_length = TG.sample(model, TG.Deterministic(), seeds; nt = 1000);
nothing
```

```@example FCUP
f, sc = @time TG.plot_odf(model; n_sphere = 500, radius = 0.3, st = 2);
TG.plot_streamlines!(sc, streamlines[1:3, 1:1:end, 1:10:Nmc])
f
```