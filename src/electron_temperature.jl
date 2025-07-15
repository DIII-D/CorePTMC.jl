## -------- PLASMA BACKGROUND ELECTRON TEMPERATURE MODELS -------- ## 

# CONSTANT Tₑ
struct TₑSurfaceModel{T} <: AbstractPlasmaModel  # AbstractPlasmaModel defined in 'plasma.jl'
    Tₑ::T# Tₑ at sheath entrance 
end
(model::TₑSurfaceModel)(Ns::Int64) = TₑSurfaceModel(zeros(Float64, Ns) .+ model.Tₑ)
TₑSurfaceModel() = TₑSurfaceModel(0.0)
struct ConstantElecTempModel{T} <: AbstractPlasmaModel  # AbstractPlasmaModel defined in 'plasma.jl'
    Tₑ::T# Tₑ at sheath entrance 
end
ConstantElecTempModel() = ConstantElecTempModel(0.0)
(model::ConstantElecTempModel)(Ns::Int64) = ConstantElecTempModel(zeros(Float64, Ns) .+ model.Tₑ)
constant_elec_temp(args...; kw...) = ConstantElecTempModel(args...; kw...)
function Tₑ_model(model::Symbol) 
    model == :surface && return TₑSurfaceModel()
    error("Unknown electron temperature model: $model")
end 

function setup_model!(model::ConstantElecTempModel, i::Int)
    nothing
end

function setup_model!(model::TₑSurfaceModel, Tₑ::Float64, i::Int)
    model.Tₑ[i] = Tₑ
    nothing
end

Tₑ_func(p::ParticlePosition, Tₑ_model::ConstantElecTempModel, i::Int64) = Tₑ_model.Tₑ 

#function update_Tₑ!(Tₑ::ElectronTemperature, p::ParticlePosition, Tₑ_model::ConstantElecTempModel, i::Int64)
#    Tₑ.value[i] = Tₑ_model.Tₑ
#end

#function update_Tₑ!(Tₑ::ElectronTemperature, p, Tₑ_model::ConstantElecTempModel)
#    Tₑ.value .= Tₑ_model.Tₑ
#end

# LINEAR Tₑ
struct Linear1DElecTempModel <: AbstractPlasmaModel
    Tₑ::Float64           # Tₑ at sheath entrance 
    T_wall::Float64       # Tₑ at wall 
    δ_sheath::Float64     # sheath width
    z_surface::Float64    # surface coordinate
end

# Constructor wrapper
linear1D_temp(args...; kw...) = Linear1DElecTempModel(args...; kw...)

# heaviside function defined in 'utils.jl'
Tₑ_func(p::ParticlePosition, Tₑ_model::Linear1DElecTempModel, i::Int64) = Tₑ_model.T_wall + (Tₑ_model.Tₑ - Tₑ_model.T_wall) / Tₑ_model.δ_sheath * p.z[i] - heaviside(p.z[i] - Tₑ_model.δ_sheath) * 
(Tₑ_model.Tₑ - Tₑ_model.T_wall) / Tₑ_model.δ_sheath * ( p.z[i] - Tₑ_model.δ_sheath )

Tₑ_func(p::ParticlePosition, Tₑ_model::Linear1DElecTempModel) = @. Tₑ_model.T_wall + (Tₑ_model.Tₑ - Tₑ_model.T_wall) / Tₑ_model.δ_sheath * p.z - heaviside(p.z - Tₑ_model.δ_sheath) * 
(Tₑ_model.Tₑ - Tₑ_model.T_wall) / Tₑ_model.δ_sheath * ( p.z[i] - Tₑ_model.δ_sheath )

function update_Tₑ!(Tₑ::ElectronTemperature,  p::ParticlePosition, Tₑ_model::AbstractPlasmaModel, i::Int64)
    Tₑ.value[i] = Tₑ_func(p, Tₑ_model, i)
end

function update_Tₑ!(Tₑ::ElectronTemperature,  p::ParticlePosition, Tₑ_model::AbstractPlasmaModel)
    @. Tₑ.value = Tₑ_func(p, Tₑ_model)
end

compute_Tₑ(model::TₑSurfaceModel, Te, i::Int64) = model.Tₑ[i]