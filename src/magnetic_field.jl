abstract type AbstractMagneticFieldModel end

struct ConstantMagneticField <: AbstractMagneticFieldModel
    B₀::Float64 # Magnetic field magnitude
    α::Float64 # Pitch angle of magnetic field
    B::Vector3D{Float64,MagneticFieldType}
end
B₀(b::MagneticField) = vector_norm(b)
inclined_magnetic_field(B₀, α) = MagneticField(B₀ * cos(α), 0.0, -B₀ * sin(α))
Bx(p::ParticlePosition, model::ConstantMagneticField) = model.B.x
By(p::ParticlePosition, model::ConstantMagneticField, args...) = model.B.y
Bz(p::ParticlePosition, model::ConstantMagneticField) = model.B.z

Bx(p::ParticlePosition, B::MagneticField, i) = B.x
By(p::ParticlePosition, B::MagneticField, i) = B.y
Bz(p::ParticlePosition, B::MagneticField, i) = B.z

function update_B!(B::MagneticField, p::ParticlePosition, model::AbstractMagneticFieldModel)
    @. B.x = Bx(p, model)
    @. B.y = By(p, model)
    @. B.z = Bz(p, model)  # Uniform magnetic field along z-axis
end

function update_B!(B::MagneticField, p::ParticlePosition, model::Union{MagneticField,AbstractMagneticFieldModel}, i)
    B.x[i] = Bx(p, model, i)
    B.y[i] = By(p, model, i)
    B.z[i] = Bz(p, model, i)  # Uniform magnetic field along z-axis
end