# Simple example of tractography

This is a more advanced tutorial because we want to show how to apply a mask.

```@example BASIC
import Tractography as TG

# define the model for TMC
model = TG.TMC(Δt = 0.125f0,
            odfdata = TG.ODFData((@__DIR__) * "/../../examples/cross-fod.nii.gz"),
            cone = TG.Cone(15f0),
            proba_min = 0.005f0,
            )
```
