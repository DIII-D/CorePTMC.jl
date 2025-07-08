## -------- PLASMA BACKGROUND MAIN ION TEMPERATURE MODELS -------- ## 

# CONSTANT Tᵢ
struct ConstantMainIonTempModel <: AbstractPlasmaModel  # AbstractPlasmaModel defined in 'plasma.jl'
    Tᵢ::Float64 # Tᵢ at sheath entrance 
end
constant_mainion_temp(args...; kw...) = ConstantMainIonTempModel(args...; kw...)

function update_Tᵢ!(Tᵢ::MainIonTemperature, p, Tᵢ_model::ConstantMainIonTempModel, i)
    Tᵢ.value[i] = Tᵢ_model.Tᵢ
end