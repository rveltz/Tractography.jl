using Tractography, Test
const TG = Tractography


TG.∂θro_sh(0.,0.,1,1)
TG.∂ϕro_sh(0.,0.,1,1)
TG.get_vector_of_sh([(0.,0.), (0.1, 0.1)], 8, 0)
TG.get_vector_of_sh([(0.,0.), (0.1, 0.1)], 8, 1)
TG.get_vector_of_sh([(0.,0.), (0.1, 0.1)], 8, 2)

# test of the value of real spherical harmonics
# https://userdocs.mrtrix.org/en/latest/concepts/spherical_harmonics.html
# https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
# https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
# careful, $Y_l^m \neq Y_{l,m}$
my_ylm(t,p,ii) = TG.get_vector_of_sh([(t,p)],8)[ii]

θ = 0.1f0; ϕ = 0.2f0
@test my_ylm(θ, ϕ, 1) ≈ sqrt(1/4/pi)

# recall that real spherical harmonics Y2,-2 ≡ sqrt(15/pi)/4*sin(2ϕ)*sin(θ)^2
fact_sh(l,m) = sqrt((2l+1)//4/pi * factorial(l-m)//factorial(l+m) )
@test fact_sh(2,-2) * 3/24 * sqrt(2) ≈ sqrt(15/pi)/4

@test my_ylm(θ, ϕ, 2) ≈ sqrt(15/pi)/4*sin(2ϕ)*sin(θ)^2
@test my_ylm(θ, ϕ, 2) ≈ 0.54627424*  sin(2ϕ)*(sin(θ)^2)

@test my_ylm(θ, ϕ, 3) ≈ -0.54627424*sin(2θ)*sin(ϕ)
@test my_ylm(θ, ϕ, 3) ≈ -sqrt(15/pi)/4*sin(2θ)*sin(ϕ)

@test my_ylm(θ, ϕ, 4) ≈ -0.31539157 + 0.94617474(cos(θ)^2)

@test my_ylm(θ, ϕ, 5) ≈ -sqrt(15/pi)/4*cos(ϕ)*sin(2θ)
@test my_ylm(θ, ϕ, 6) ≈ 0.54627424(sin(θ)^2)*cos(2ϕ)

@test my_ylm(θ, ϕ, 7) ≈ 0.6258357sin(4ϕ)*(sin(θ)^4)
@test my_ylm(θ, ϕ, 8) ≈ -1.7701306cos(θ)*sin(3ϕ)*(sin(θ)^3)
@test my_ylm(θ, ϕ, 9) ≈ -0.47308737sin(2ϕ)*(sin(θ)^2) + 3.3116114(cos(θ)^2)*sin(2ϕ)*(sin(θ)^2)
@test my_ylm(θ, ϕ, 10) ≈ 1.0035698sin(2θ)*sin(ϕ) - 4.683326(cos(θ)^3)*sin(θ)*sin(ϕ)
@test my_ylm(θ, ϕ, 11) ≈ 0.31735665 - 3.1735663(cos(θ)^2) + 3.7024941(cos(θ)^4)
@test my_ylm(θ, ϕ, 12) ≈ 1.0035698cos(ϕ)*sin(2θ) - 4.683326(cos(θ)^3)*cos(ϕ)*sin(θ)
@test my_ylm(θ, ϕ, 13) ≈ -0.47308737(sin(θ)^2)*cos(2ϕ) + 3.3116114(cos(θ)^2)*(sin(θ)^2)*cos(2ϕ)
@test my_ylm(θ, ϕ, 14) ≈ -1.7701306cos(θ)*(sin(θ)^3)*cos(3ϕ)
@test my_ylm(θ, ϕ, 15) ≈ 0.6258357(sin(θ)^4)*cos(4ϕ)