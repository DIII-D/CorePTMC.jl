## -------- PLASMA BACKGROUND ELECTRON DENSITY MODELS -------- ##

struct BoltzmannElectronDensityModel{T} <: AbstractPlasmaModel
    nₑ₀::T
    Tₑ::T
end
BoltzmannElectronDensityModel() = BoltzmannElectronDensityModel(0.0, 0.0)
(model::BoltzmannElectronDensityModel)(Ns::Int64) = BoltzmannElectronDensityModel(zeros(Float64, Ns) .+ model.nₑ₀, zeros(Float64, Ns) .+ model.Tₑ)

function setup_model!(model::BoltzmannElectronDensityModel,  Tₑ::Float64, nₑ₀, i::Int)
    model.nₑ₀[i], model.Tₑ[i] = nₑ₀, Tₑ
    nothing
end

struct ConstantElectronDensityModel <: AbstractPlasmaModel
    nₑ₀::Float64
    #Tₑ::Float64
end

constant_elec_density(args...; kw...) = ConstantElectronDensityModel(args...; kw...)
function update_nₑ!(nₑ::ElectronDensity, p, ϕ::ElectricPotential, nₑ_model::ConstantElectronDensityModel, i)
    nₑ.value[i] = nₑ_model.nₑ₀
end

function update_nₑ!(nₑ::ElectronDensity, p, ϕ::ElectricPotential, nₑ_model::ConstantElectronDensityModel)
    nₑ.value .= nₑ_model.nₑ₀
end
nₑ(model::BoltzmannElectronDensityModel, ϕ::Float64, i::Int64) = model.nₑ₀[i] * exp(ϕ / model.Tₑ[i])
boltzmann_elec_density(args...; kw...) = BoltzmannElectronDensityModel(args...; kw...)
function update_nₑ!(nₑ::ElectronDensity, ϕ::ElectricPotential, nₑ_model::BoltzmannElectronDensityModel)
    nₑ.value = nₑ_boltzmann(ϕ, nₑ_model.nₑ₀, nₑ_model.Tₑ)
end

function update_nₑ!(nₑ::ElectronDensity, p, ϕ::ElectricPotential, nₑ_model::BoltzmannElectronDensityModel, i)
    nₑ.value[i] = nₑ_boltzmann(ϕ[i], nₑ_model.nₑ₀, nₑ_model.Tₑ)
end
compute_nₑ(args...; kw...) = nₑ(args...; kw...)
nₑ_boltzmann(ϕ, nₑ₀, Tₑ) = @. nₑ₀ * exp(ϕ / Tₑ)
nₑ_model(model::Symbol) = model == :boltzmann ? BoltzmannElectronDensityModel() : error("Unknown electron density model: $model")
#TODO --- Luca: we could add a kinetic tail to the electron distribution and have f_e instead of n_e but this is for mid-term/long-term future development ---  