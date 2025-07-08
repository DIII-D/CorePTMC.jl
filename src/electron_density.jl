## -------- PLASMA BACKGROUND ELECTRON DENSITY MODELS -------- ##

struct BoltzmannElectronDensityModel <: AbstractPlasmaModel
    nₑ₀::Float64
    Tₑ::Float64
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

boltzmann_elec_density(args...; kw...) = BoltzmannElectronDensityModel(args...; kw...)
function update_nₑ!(nₑ::ElectronDensity, ϕ::ElectricPotential, nₑ_model::BoltzmannElectronDensityModel)
    nₑ.value = nₑ_boltzmann(ϕ, nₑ_model.nₑ₀, nₑ_model.Tₑ)
end

function update_nₑ!(nₑ::ElectronDensity, p, ϕ::ElectricPotential, nₑ_model::BoltzmannElectronDensityModel, i)
    nₑ.value[i] = nₑ_boltzmann(ϕ[i], nₑ_model.nₑ₀, nₑ_model.Tₑ)
end

nₑ_boltzmann(ϕ, nₑ₀, Tₑ) = @. nₑ₀ * exp(ϕ / Tₑ)

#TODO --- Luca: we could add a kinetic tail to the electron distribution and have f_e instead of n_e but this is for mid-term/long-term future development ---  