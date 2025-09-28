using Tractography, StaticArrays
const TG = Tractography

_ni = ones(145, 175, 145, 45);
S = 1.0 * I(4) .+ 0.1randn(4,4) |> Array
Si = pinv(S)
T = zeros(3)
model = TG.TMC(Δt = 0.1, proba_min = 0.0,
            odfdata = ODFData("bla.txt",
                                     _ni,
                                     S,
                                     T
                                     ),
            )

x = SVector(rand(3)...)
y = S * SVector(x..., 1)
@test TG.transform(model.odfdata, x) == y
@test TG.transform_inv(model.odfdata, y)[1:3] ≈ x

TG.get_range(model.odfdata)