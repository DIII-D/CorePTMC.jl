## -------- PLASMA BACKGROUND MAIN ION DENSITY MODELS -------- ##

struct BoltzmannMainIonDensityModel <: AbstractPlasmaModel
    nᵢ₀::Float64
    Tᵢ::Float64
end

struct ConstantMainIonDensityModel <: AbstractPlasmaModel
    nᵢ₀::Float64
end

constant_mainion_density(args...; kw...) = ConstantMainIonDensityModel(args...; kw...)
function update_nᵢ!(nᵢ::MainIonDensity, p, ϕ::ElectricPotential, nᵢ_model::ConstantMainIonDensityModel, i)
    nᵢ.value[i] = nᵢ_model.nᵢ₀
end

function update_nᵢ!(nᵢ::MainIonDensity, p, ϕ::ElectricPotential, nᵢ_model::ConstantMainIonDensityModel)
    nᵢ.value .= nᵢ_model.nᵢ₀
end

boltzmann_mainion_density(args...; kw...) = BoltzmannMainIonDensityModel(args...; kw...)
function update_nᵢ!(nᵢ::MainIonDensity, ϕ::ElectricPotential, nᵢ_model::BoltzmannMainIonDensityModel)
    nᵢ.value = nᵢ_boltzmann(ϕ, nᵢ_model.nᵢ₀, nᵢ_model.Tᵢ)
end

function update_nₑ!(nᵢ::MainIonDensity, p, ϕ::ElectricPotential, nᵢ_model::BoltzmannMainIonDensityModel, i)
    nᵢ.value[i] = nᵢ_boltzmann(ϕ[i], nᵢ_model.nᵢ₀, nᵢ_model.Tₑ)
end

nᵢ_boltzmann(ϕ, nᵢ₀, Tₑ) = @. nᵢ₀ * exp(ϕ / Tₑ)