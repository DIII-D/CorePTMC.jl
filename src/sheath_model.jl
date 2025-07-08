abstract type AbstractSheathModel end
abstract type Abstract0DSheathModel <: AbstractSheathModel end
abstract type Abstract1DSheathModel <: AbstractSheathModel end


# NO SHEATH MODEL

struct NoSheathModel <: Abstract0DSheathModel
    ϕ_0::Float64         # plasma potential
    ϕ_wall::Float64      # wall potential
    δ_sheath::Float64    # sheath width (assumed proportional to ρᵢ: δ_sheath = f_sheath * ρᵢ)
    z_surface::Float64
end

no_sheath_model() = NoSheathModel(0.0, 0.0, 0.0, 0.0) #renaming and creation of model parameters

# No electric field if potential is constant
Ex(p::ParticlePosition, model::Abstract0DSheathModel, args...) = 0.0
Ey(p::ParticlePosition, model::Abstract0DSheathModel, args...) = 0.0
Ez(p::ParticlePosition, model::Abstract0DSheathModel, args...) = 0.0

ϕ(p::ParticlePosition, model::NoSheathModel, i) = model.ϕ_wall


# If sheath model is 1D - regardless of particle position, sheath model type and any other argument -> Ex = Ey = 0
Ex(p::ParticlePosition, model::Abstract1DSheathModel, args...) = 0.0
Ey(p::ParticlePosition, model::Abstract1DSheathModel, args...) = 0.0

# potential drop with respect to reference (plasma potential)
# Δϕ = ϕ_wall - ϕ_0
# Up to date, ϕ_0 ≡ 0. Future developments will allow the user to modify the plasma potential 

# LINEAR 1D POTENTIAL DROP

struct Linear1DSheathModel <: Abstract1DSheathModel
    ϕ_0::Float64 # plasma potential
    ϕ_wall::Float64 # wall potential
    δ_sheath::Float64
    z_surface::Float64
end

linear_sheath_model(Tₑ, ρᵢ; Λ=-3, f_sheath=5.0) = Linear1DSheathModel(0.0, Λ * Tₑ, f_sheath * ρᵢ, 0.0)


function ϕ(p::ParticlePosition, model::Linear1DSheathModel, i)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z[i] - model.z_surface

    return Δϕ * (1 - Δz / model.δ_sheath) * heaviside(model.δ_sheath - Δz)
end

function ϕ(p::ParticlePosition, model::Linear1DSheathModel)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z .- model.z_surface

    return @. Δϕ * (1 - Δz / model.δ_sheath) * heaviside(model.δ_sheath - Δz)
end

# calculating electric potential at particle's position 
function Ez(p::ParticlePosition, model::Linear1DSheathModel, i)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z[i] - model.z_surface

    return Δϕ / model.δ_sheath * heaviside(model.δ_sheath - Δz)
end

function Ez(p::ParticlePosition, model::Linear1DSheathModel)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z .- model.z_surface

    return @. Δϕ / model.δ_sheath * heaviside(model.δ_sheath - Δz)
end


## QUADRATIC 1D POTENTIAL DROP 
struct Quadratic1DSheathModel <: Abstract1DSheathModel
    ϕ_0::Float64
    ϕ_wall::Float64
    δ_sheath::Float64
    z_surface::Float64
end

quadratic_sheath_model(Tₑ, ρᵢ; Λ=-3, f_sheath=5.0) = Quadratic1DSheathModel(0.0, Λ * Tₑ, f_sheath * ρᵢ, 0.0)

function ϕ(p::ParticlePosition, model::Quadratic1DSheathModel, i)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z[i] - model.z_surface

    return Δϕ * (Δz / model.δ_sheath - 1)^2 * heaviside(model.δ_sheath - Δz)
end

function ϕ(p::ParticlePosition, model::Quadratic1DSheathModel)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z .- model.z_surface

    return @. Δϕ * (Δz / model.δ_sheath - 1)^2 * heaviside(model.δ_sheath - Δz)
end

function Ez(p::ParticlePosition, model::Quadratic1DSheathModel, i)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z[i] - model.z_surface

    return 2 * Δϕ / model.δ_sheath * (1 - Δz / model.δ_sheath) * heaviside(model.δ_sheath - Δz)
end

function Ez(p::ParticlePosition, model::Quadratic1DSheathModel)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z .- model.z_surface

    return @. 2 * Δϕ / model.δ_sheath * (1 - Δz / model.δ_sheath) * heaviside(model.δ_sheath - Δz)
end


# MISCELLANEOUS FUNCTIONS
heaviside(x) = @. (1.0 + sign(x)) / 2


## ---- exponential sheath model ---- # from J. Guterl NME 2021 W redeposition
struct Exponential1DSheathModel <: Abstract1DSheathModel
    ϕ_0::Float64
    ϕ_wall::Float64
    δ_sheath::Float64
    z_surface::Float64
    β_sheath::Float64
end


exponential_sheath_model(Tₑ, ρᵢ; Λ=-3, f_sheath=5.0, β_sheath=4.0) = Exponential1DSheathModel(0.0, Λ * Tₑ, f_sheath * ρᵢ, 0.0, β_sheath)

function ϕ(p::ParticlePosition, model::Exponential1DSheathModel, i)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z[i] - model.z_surface

    return Δϕ * exp(-model.β_sheath * Δz / model.δ_sheath)
end

function ϕ(p::ParticlePosition, model::Exponential1DSheathModel)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z .- model.z_surface

    return @. Δϕ * exp(-model.β_sheath * Δz / model.δ_sheath)
end

function Ez(p::ParticlePosition, model::Exponential1DSheathModel, i)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z[i] - model.z_surface

    return Δϕ * model.β_sheath / model.δ_sheath * exp(-model.β_sheath * Δz / model.δ_sheath)
end

function Ez(p::ParticlePosition, model::Exponential1DSheathModel)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z .- model.z_surface

    return @. Δϕ * model.β_sheath / model.δ_sheath * exp(-model.β_sheath * Δz / model.δ_sheath)
end

function grad_Ez(p::ParticlePosition, model::Exponential1DSheathModel)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z .- model.z_surface

    return @. - Δϕ * (model.β_sheath / model.δ_sheath) ^ 2 * exp(-model.β_sheath * Δz / model.δ_sheath)
end

# ---- Borodkina sheath model ---- #

# paper1) Borodkina et al. 2015 Russian Physics Journal DOI 10.1007/s11182-015-0518-5
# paper2) Borodkina et al. 2016 Contrib. Plasma Physics DOI 10.1002/ctpp.201610032

struct BorodkinaSheathModel <: Abstract1DSheathModel
    ϕ_0::Float64
    ϕ_wall::Float64
    Tₑ::Float64
    α::Float64
    ρᵢ::Float64
    K::Float64
    z_surface::Float64
    r_d::Float64        # Debye radius
end

borodkina_sheath_model(Tₑ, ρᵢ, α, r_d; Λ=-3, bias=0.0, K=2.0) = BorodkinaSheathModel(0.0, (Λ + bias) * Tₑ, Tₑ, α , ρᵢ, K, 0.0, r_d)

function ϕ(p::ParticlePosition, model::BorodkinaSheathModel, i)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z[i] - model.z_surface
    r_d = model.r_d
    Tₑ = model.Tₑ
    ρᵢ = model.ρᵢ
    ξ = Δz / r_d
    
    Λ_w = model.ϕ_wall / Tₑ        # total potential drop including bias

    α_B = pi / 2 - model.α               # Borodkina's definition of magnetic field poloidal angle 
    α_crit = acos( exp(Λ_w) )      # critical angle below which Debye sheath disappears

    if α_B >= α_crit
        DS = false
    else
        DS = true
    end

    if DS
        Λ_mps = log( cos( α_B ) )   # Paper1 - eq. 6 (potential drop in mps)
        dn_mps = - Λ_mps / ( ρᵢ / r_d * sin( α_B ) ) ^ 2     # Paper1 - eq. 12 (difference between ion and electron densities at mps entrance)
        C_1 = -dn_mps * Λ_mps - 6 * cos( α_B )               # Paper1 - eq. 20 (Integration constant at mps entrance)

        # potential drop curve parameters (a, Q): Paper1 - set of eqs. 23 factor is a variable written to simplify coding of (a, Q)
        factor = sqrt( 2.0 * exp( Λ_w ) + 4.0 * cos( α_B ) * sqrt( 1.0 - ( Λ_w - Λ_mps ) ) + C_1 )
        a = ( sqrt( -dn_mps * Λ_mps ) - factor ) ./ ( Λ_w - Λ_mps )
        Q = factor / a
        ξ_mps = -1 / a * log( ( Λ_w - Λ_mps + Q ) / Q )       # Paper1 - eq. 24 (normalized Debye Sheath width)
        L_mps = model.K * ρᵢ  / r_d * sin( α_B )              # Paper2 - eq. 2  (length of mps)
        ϕ_ds = Δϕ + Q * Tₑ * ( 1 - exp( -a * ξ ) )            # Paper1 - eq. 22 (debye sheath potential drop)
        ϕ_mps = Λ_mps * Tₑ * exp( -2 * ( ξ - ξ_mps ) / L_mps) # Paper2 - eq. 2  (mps sheath potential drop)

        if ξ < ξ_mps
            return ϕ_ds
        else
            return ϕ_mps
        end

    else
        Λ_mps = Λ_w
        ξ_mps = 0.0
        L_mps = model.K * ρᵢ  / r_d * sin( α_B )              # Paper2 - eq. 2 (length of mps)
        ϕ_mps = Λ_mps * Tₑ * exp( -2 * ( ξ - ξ_mps ) / L_mps) # Paper2 - eq. 2 (mps sheath potential drop)

        return ϕ_mps
    end

end

function ϕ(p::ParticlePosition, model::BorodkinaSheathModel)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z .- model.z_surface
    r_d = model.r_d
    Tₑ = model.Tₑ
    ρᵢ = model.ρᵢ
    ξ = @. Δz / r_d
    
    ϕ = zeros(length(p.z))         # initalizing potential drop

    Λ_w = model.ϕ_wall / Tₑ        # total potential drop including bias

    α_B = pi / 2 - model.α         # Borodkina's definition of magnetic field poloidal angle 
    α_crit = acos( exp(Λ_w) )      # critical angle below which Debye sheath disappears

    if α_B >= α_crit
        DS = false
    else
        DS = true
    end

    if DS
        Λ_mps = log( cos( α_B ) )   # Paper1 - eq. 6 (potential drop in mps)
        dn_mps = - Λ_mps / ( ρᵢ / r_d * sin( α_B ) ) ^ 2     # Paper1 - eq. 12 (difference between ion and electron densities at mps entrance)
        C_1 = -dn_mps * Λ_mps - 6 * cos( α_B )               # Paper1 - eq. 20 (Integration constant at mps entrance)

        # potential drop curve parameters (a, Q): Paper1 - set of eqs. 23 factor is a variable written to simplify coding of (a, Q)
        factor = sqrt( 2.0 * exp( Λ_w ) + 4.0 * cos( α_B ) * sqrt( 1.0 - ( Λ_w - Λ_mps ) ) + C_1 )
        a = ( sqrt( -dn_mps * Λ_mps ) - factor ) ./ ( Λ_w - Λ_mps )
        Q = factor / a
        ξ_mps = -1 / a * log( ( Λ_w - Λ_mps + Q ) / Q )       # Paper1 - eq. 24 (normalized Debye Sheath width)
        L_mps = model.K * ρᵢ  / r_d * sin( α_B )              # Paper2 - eq. 2  (length of mps)
        
        ϕ[ξ .< ξ_mps] = @. Δϕ + Q * Tₑ * ( 1 - exp( -a * ξ[ξ .< ξ_mps] ) )           # Paper1 - eq. 22 (debye sheath potential drop)
        ϕ[ξ .> ξ_mps] = @. Λ_mps * Tₑ * exp( -2 * ( ξ[ξ .> ξ_mps] - ξ_mps ) / L_mps) # Paper2 - eq. 2 (mps sheath potential drop)

    else
        Λ_mps = Λ_w
        ξ_mps = 0.0
        L_mps = model.K * ρᵢ  / r_d * sin( α_B )              # Paper2 - eq. 2 (length of mps)
        ϕ = @. Λ_mps * Tₑ * exp( -2 * ( ξ - ξ_mps ) / L_mps)  # Paper2 - eq. 2 (mps sheath potential drop)
    end

    return ϕ
end



function Ez(p::ParticlePosition, model::BorodkinaSheathModel, i)
    
    Δz = p.z[i] - model.z_surface
    r_d = model.r_d
    Tₑ = model.Tₑ
    ρᵢ = model.ρᵢ
    ξ = Δz / r_d
    
    Λ_w = model.ϕ_wall / Tₑ        # total potential drop including bias

    α_B = pi / 2 - model.α         # Borodkina's definition of magnetic field poloidal angle 
    α_crit = acos( exp(Λ_w) )      # critical angle below which Debye sheath disappears

    if α_B >= α_crit
        DS = false
    else
        DS = true
    end

    if DS
        Λ_mps = log( cos( α_B ) )   # Paper1 - eq. 6 (potential drop in mps)
        dn_mps = - Λ_mps / ( ρᵢ / r_d * sin( α_B ) ) ^ 2     # Paper1 - eq. 12 (difference between ion and electron densities at mps entrance)
        C_1 = -dn_mps * Λ_mps - 6 * cos( α_B )               # Paper1 - eq. 20 (Integration constant at mps entrance)

        # potential drop curve parameters (a, Q): Paper1 - set of eqs. 23 factor is a variable written to simplify coding of (a, Q)
        factor = sqrt( 2.0 * exp( Λ_w ) + 4.0 * cos( α_B ) * sqrt( 1.0 - ( Λ_w - Λ_mps ) ) + C_1 )
        a = ( sqrt( -dn_mps * Λ_mps ) - factor ) ./ ( Λ_w - Λ_mps )
        Q = factor / a
        ξ_mps = -1 / a * log( ( Λ_w - Λ_mps + Q ) / Q )                         # Paper1 - eq. 24 (normalized Debye Sheath width)
        L_mps = model.K * ρᵢ  / r_d * sin( α_B )                                # Paper2 - eq. 2  (length of mps)
        E_ds = -a * Q * Tₑ / r_d * exp( -a * ξ )                                # derivative of Paper1 - eq. 22 (debye sheath electric field)
        E_mps = Λ_mps * Tₑ * 2 / L_mps / r_d * exp(- 2 / L_mps * ( ξ - ξ_mps )) # derivative of Paper2 - eq. 2 (mps sheath electric field)

        if ξ < ξ_mps
            return E_ds
        else
            return E_mps
        end

    else
        Λ_mps = Λ_w
        ξ_mps = 0.0
        L_mps = model.K * ρᵢ  / r_d * sin( α_B )                                # Paper2 - eq. 2 (length of mps)
        E_mps = Λ_mps * Tₑ * 2 / L_mps / r_d * exp(- 2 / L_mps * ( ξ - ξ_mps )) # derivative of Paper2 - eq. 2 (mps sheath electric field)

        return E_mps
    end

end


function Ez(p::ParticlePosition, model::BorodkinaSheathModel)
    Δz = p.z .- model.z_surface
    r_d = model.r_d
    Tₑ = model.Tₑ
    ρᵢ = model.ρᵢ
    ξ = @. Δz / r_d
    
    E = zeros(length(p.z))         # initalizing potential drop

    Λ_w = model.ϕ_wall / Tₑ        # total potential drop including bias

    α_B = pi / 2 - model.α         # Borodkina's definition of magnetic field poloidal angle 
    α_crit = acos( exp(Λ_w) )      # critical angle below which Debye sheath disappears

    if α_B >= α_crit
        DS = false
    else
        DS = true
    end

    if DS
        Λ_mps = log( cos( α_B ) )                            # Paper1 - eq. 6 (potential drop in mps)
        dn_mps = - Λ_mps / ( ρᵢ / r_d * sin( α_B ) ) ^ 2     # Paper1 - eq. 12 (difference between ion and electron densities at mps entrance)
        C_1 = -dn_mps * Λ_mps - 6 * cos( α_B )               # Paper1 - eq. 20 (Integration constant at mps entrance)

        # potential drop curve parameters (a, Q): Paper1 - set of eqs. 23 factor is a variable written to simplify coding of (a, Q)
        factor = sqrt( 2.0 * exp( Λ_w ) + 4.0 * cos( α_B ) * sqrt( 1.0 - ( Λ_w - Λ_mps ) ) + C_1 )
        a = ( sqrt( -dn_mps * Λ_mps ) - factor ) ./ ( Λ_w - Λ_mps )
        Q = factor / a
        ξ_mps = -1 / a * log( ( Λ_w - Λ_mps + Q ) / Q )       # Paper1 - eq. 24 (normalized Debye Sheath width)
        L_mps = model.K * ρᵢ  / r_d * sin( α_B )              # Paper2 - eq. 2  (length of mps)
        
        E[ξ .< ξ_mps] = @. -a * Q * Tₑ / r_d * exp( -a * ξ[ξ .< ξ_mps] )                               # derivative of Paper1 - eq. 22 (debye sheath electric field)
        E[ξ .> ξ_mps] = @. Λ_mps * Tₑ * 2 / L_mps / r_d * exp(- 2 / L_mps * ( ξ[ξ .> ξ_mps] - ξ_mps )) # derivative of Paper2 - eq. 2 (mps sheath electric field)

    else
        Λ_mps = Λ_w
        ξ_mps = 0.0
        L_mps = model.K * ρᵢ  / r_d * sin( α_B )              # Paper2 - eq. 2 (length of mps)
        E = @. Λ_mps * Tₑ * 2 / L_mps / r_d * exp(- 2 / L_mps * ( ξ[ξ .> ξ_mps] - ξ_mps )) # derivative of Paper2 - eq. 2 (mps sheath electric field)
    end

    return E
end

# ----- Brooks' model ----- #

# N. Brooks Modeling of sputtering erosion/redeposition—status and implications 
# for fusion design Fusion Engineering and Design 60 (2002) 515–526 


struct Brooks1DSheathModel <: Abstract1DSheathModel
    ϕ_0::Float64
    ϕ_wall::Float64
    α::Float64
    δ_sheath::Float64
    z_surface::Float64
    r_d::Float64        # Debye radius
    Tₑ::Float64
end

brooks_sheath_model(Tₑ, ρᵢ, α, r_d; Λ=-3, f_sheath=1.0, bias = 0.0) = Brooks1DSheathModel(0.0, Λ * Tₑ + bias, α, f_sheath * ρᵢ, 0.0, r_d, Tₑ)

function ϕ(p::ParticlePosition, model::Brooks1DSheathModel, i)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z[i] - model.z_surface
    
    α_B = π/2 - model.α   # Brooks' definition of magnetic field poloidal angle 

    Δϕ₂ = model.Tₑ * log(cos( α_B )) # fraction of Δϕ in Chodura sheath. f_d(α_B = 0) = 0 -> only Debye sheath. f_d(α_B = π/2) = 
    
    if abs(Δϕ) < abs(Δϕ₂) # if α_B is too shallow (close to 90 deg) all potential drop happens inside the chodura sheath
        Δϕ₁ = 0  # no debye sheath
        Δϕ₂ = Δϕ 
    else
        Δϕ₁ = Δϕ - Δϕ₂
    end

    f_d = abs(Δϕ₁ / Δϕ)  # fraction of potential drop in Debye sheath

    return Δϕ * (f_d * exp(-Δz / (2 * model.r_d)) + (1 - f_d) * exp(-Δz / model.δ_sheath) )
end

function ϕ(p::ParticlePosition, model::Brooks1DSheathModel)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z .- model.z_surface
    
    α_B = π/2 - model.α   # Brooks' definition of magnetic field poloidal angle 

    Δϕ₂ = model.Tₑ * log(cos( α_B )) # fraction of Δϕ in Chodura sheath. f_d(α_B = 0) = 0 -> only Debye sheath. f_d(α_B = π/2) = 
    
    #@show Δϕ₂

    #@show abs(Δϕ), abs(Δϕ₂)

    if abs(Δϕ) < abs(Δϕ₂) # if α_B is too shallow (close to 90 deg) all potential drop happens inside the chodura sheath
        Δϕ₁ = 0  # no debye sheath
        Δϕ₂ = Δϕ 
    else
        Δϕ₁ = Δϕ - Δϕ₂
    end

    f_d = abs(Δϕ₁ / Δϕ)  # fraction of potential drop in Debye sheath

    #@show f_d

    return @. Δϕ * (f_d * exp(-Δz / (2 * model.r_d)) + (1 - f_d) * exp(-Δz / model.δ_sheath) )
end

function Ez(p::ParticlePosition, model::Brooks1DSheathModel, i)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z[i] - model.z_surface

    α_B = π/2 - model.α   # Brooks' definition of magnetic field poloidal angle 

    Δϕ₂ = model.Tₑ * log(cos( α_B )) # fraction of Δϕ in Chodura sheath. f_d(α_B = 0) = 0 -> only Debye sheath. f_d(α_B = π/2) = 
    
    if abs(Δϕ) < abs(Δϕ₂) # if α_B is too shallow (close to 90 deg) all potential drop happens inside the chodura sheath
        Δϕ₁ = 0  # no debye sheath
        Δϕ₂ = Δϕ 
    else
        Δϕ₁ = Δϕ - Δϕ₂
    end

    f_d = abs(Δϕ₁ / Δϕ)  # fraction of potential drop in Debye sheath

    return Δϕ * (f_d / (2 * model.r_d) * exp(-Δz / (2 * model.r_d)) + (1 - f_d) / model.δ_sheath * exp(-Δz / model.δ_sheath) )
end

function Ez(p::ParticlePosition, model::Brooks1DSheathModel)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z .- model.z_surface

    α_B = π/2 - model.α   # Brooks' definition of magnetic field poloidal angle 

    Δϕ₂ = model.Tₑ * log(cos( α_B )) # fraction of Δϕ in Chodura sheath. f_d(α_B = 0) = 0 -> only Debye sheath. f_d(α_B = π/2) = 
    
    if abs(Δϕ) < abs(Δϕ₂) # if α_B is too shallow (close to 90 deg) all potential drop happens inside the chodura sheath
        Δϕ₁ = 0  # no debye sheath
        Δϕ₂ = Δϕ 
    else
        Δϕ₁ = Δϕ - Δϕ₂
    end

    f_d = abs(Δϕ₁ / Δϕ)  # fraction of potential drop in Debye sheath

    return @. Δϕ * (f_d / (2 * model.r_d) * exp(-Δz / (2 * model.r_d)) + (1 - f_d) / model.δ_sheath * exp(-Δz / model.δ_sheath) )
end

function grad_Ez(p::ParticlePosition, model::Brooks1DSheathModel)
    Δϕ = model.ϕ_wall - model.ϕ_0
    Δz = p.z .- model.z_surface

    α_B = π/2 - model.α   # Brooks' definition of magnetic field poloidal angle 

    Δϕ₂ = model.Tₑ * log(cos( α_B )) # fraction of Δϕ in Chodura sheath. f_d(α_B = 0) = 0 -> only Debye sheath. f_d(α_B = π/2) = 
    
    if abs(Δϕ) < abs(Δϕ₂) # if α_B is too shallow (close to 90 deg) all potential drop happens inside the chodura sheath
        Δϕ₁ = 0  # no debye sheath
        Δϕ₂ = Δϕ 
    else
        Δϕ₁ = Δϕ - Δϕ₂
    end

    f_d = abs(Δϕ₁ / Δϕ)  # fraction of potential drop in Debye sheath

    return @. - Δϕ * (f_d / (2 * model.r_d) ^2 .* exp(-Δz / (2 * model.r_d)) + (1 - f_d)/model.δ_sheath ^ 2 * exp(-Δz / model.δ_sheath))
end

# --------------------------------------- #
# --------------------------------------- #

# updaters
function update_E!(E::ElectricField, p::ParticlePosition, model::AbstractSheathModel)
    @. E.x = Ex(p, model)
    @. E.y = Ey(p, model)
    @. E.z = Ez(p, model)
end

function update_E!(E::ElectricField, p::ParticlePosition, model::AbstractSheathModel, i::Int64)
    E.x[i] = Ex(p, model, i)
    E.y[i] = Ey(p, model, i)
    E.z[i] = Ez(p, model, i)
end

function update_ϕ!(ϕ::ElectricPotential, p::ParticlePosition, model::AbstractSheathModel)
    @. ϕ.value = ϕ(p, model)
end

function update_ϕ!(φ::ElectricPotential, p::ParticlePosition, model::AbstractSheathModel, i)
    φ.value[i] = ϕ(p, model, i)
end