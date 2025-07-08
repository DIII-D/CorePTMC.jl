## -------- PLASMA BACKGROUND MAIN ION AVERAGE VELOCITY MODELS -------- ## 

# Bohm_criterion (sound speed) 
struct ConstantDperp <: AbstractPlasmaModel  # AbstractPlasmaModel defined in 'plasma.jl'
    D_perp::Float64 
end

constant_dperp(D_perp) = ConstantDperp(D_perp) # function to call in your input file

function update_dperp!(D_perp::AnomalousDiffusion, p, dperp_model::ConstantDperp, i)
    D_perp.value[i] = dperp_model.D_perp
end