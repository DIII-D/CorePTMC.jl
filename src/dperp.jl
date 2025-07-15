## -------- PLASMA BACKGROUND MAIN ION AVERAGE VELOCITY MODELS -------- ## 

# Bohm_criterion (sound speed) 

struct ConstantDperp{T} <: AbstractPlasmaModel  # AbstractPlasmaModel defined in 'plasma.jl'
    D_perp::T
end

D_perp_model(model::Float64) = ConstantDperp(model)

(model::ConstantDperp)(Ns::Int64) = ConstantDperp(zeros(Float64, Ns) .+ model.D_perp)
constant_dperp(D_perp) = ConstantDperp(D_perp) # function to call in your input file

function update_dperp!(D_perp::AnomalousDiffusion, p, dperp_model::ConstantDperp, i)
    D_perp.value[i] = dperp_model.D_perp
end

compute_Dperp(model::ConstantDperp, Dperp, i::Int64) = model.D_perp[i]