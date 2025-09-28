# Simple example of tractography

This is a more advanced tutorial because we want to show how to apply a mask.

```@example BASIC
import Tractography as TG

# define the model for TMC
model = TG.TMC(Δt = 0.125f0,
            odfdata = TG.ODFData("/Users/rveltz/work/prog_gd/julia/dev/dev1/Tractography/examples/cross-fod.nii.gz"),
            C = TG.Cone(15),
            proba_min = 0.005f0,
            )
```
