using Tractography, StaticArrays
const TG = Tractography

model = TG.TMC(Δt = 0.1, proba_min = 0.0,
            odfdata = TG.ODFData((@__DIR__) * "/../examples/fod-FC.nii.gz")
)

show(stdout, model.odfdata; full = true)

x = SVector(rand(3)...)
y = model.odfdata.transform.S * SVector(x..., 1)
TG.transform(model.odfdata, CartesianIndex{3}(1,1,1))
@test TG.transform(model.odfdata, x) == y
@test TG.transform_inv(model.odfdata, y)[1:3] ≈ x rtol = 1e-6

TG.get_range(model.odfdata)
TG._get_array(zeros(2))
TG._my_typeof(1)

TG.from_odf(model, 10)
mask = model.odfdata.data[:,:,:,1].>0;
TG.from_mask(model, mask, 10)