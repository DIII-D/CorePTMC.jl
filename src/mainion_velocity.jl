## -------- PLASMA BACKGROUND MAIN ION AVERAGE VELOCITY MODELS -------- ## 

# Bohm_criterion (sound speed) 
struct BohmMainIonVelModel <: AbstractPlasmaModel  # AbstractPlasmaModel defined in 'plasma.jl'
    vᵢ::Vector3D{Float64,MainIonVelocityType}
end

# call 'bohm_criterion' to charge 'BohmMainIonVelModel' structure with user value. 
bohm_criterion(α_B::Float64, cₛ::Float64) = BohmMainIonVelModel(AbstractVector3D{MainIonVelocityType}(cₛ .* [cos(α_B), 0.0, -sin(α_B)]))

vᵢx(p::ParticlePosition, model::BohmMainIonVelModel, i) = model.vᵢ.x
vᵢy(p::ParticlePosition, model::BohmMainIonVelModel, i) = model.vᵢ.y
vᵢz(p::ParticlePosition, model::BohmMainIonVelModel, i) = model.vᵢ.z

# update 'vᵢ' of particle 'i', depending on particle position and chosen model 
function update_vᵢ!(vᵢ::MainIonVelocity, p::ParticlePosition, model::AbstractPlasmaModel, i)
    vᵢ.x[i] = vᵢx(p, model, i)
    vᵢ.y[i] = vᵢy(p, model, i)
    vᵢ.z[i] = vᵢz(p, model, i)  # Uniform magnetic field along z-axis
end