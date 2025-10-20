using CairoMakie, NIfTI

model = Tractography.TMC(Δt = 0.125f0,
            odfdata = Tractography.ODFData((@__DIR__) * "/../examples/fod-FC.nii.gz"),
            cone = Tractography.Cone(15),
            proba_min = 0.005f0,
            )

streamlines, tract_length = Tractography.sample(model, 
                            Tractography.Deterministic(), 
                            Tractography.from_odf(model, 10); 
                            nt = 100);

f, sc = Tractography.plot_odf(model; n_sphere = 100, radius = 0.1, st = 2);
Tractography.plot_streamlines!(sc, streamlines)
Tractography.plot_slice(model, I = 32:32)
Tractography.plot_slice(model, J = 32:32)
Tractography.plot_slice(model, K = 2:2)
