using Tractography, Test
const TG = Tractography


TG.∂θro_sh(0.,0.,1,1)
TG.∂ϕro_sh(0.,0.,1,1)
TG.get_vector_of_sh([(0.,0.), (0.1, 0.1)], 8, 0)
TG.get_vector_of_sh([(0.,0.), (0.1, 0.1)], 8, 1)
TG.get_vector_of_sh([(0.,0.), (0.1, 0.1)], 8, 2)