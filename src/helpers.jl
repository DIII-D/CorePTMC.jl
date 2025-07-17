## -------- HELPERS STRUCTURE AND CONSTRUCTOR  -------- ##

## structure for RNG
struct RandomGenerators{I,J}
    iz::I
    v₀::J
end

RandomGenerators(Np::Int64; seed=123) = RandomGenerators(MersenneTwister(seed), MersenneTwister(seed))

# structure for Helpers
struct Helpers{R<:RandomGenerators}
    random_generators::R
    α::ScalarField{Vector{Float64},GenericType}
    χ_iz::ScalarField{Vector{Float64},RandomGeneratorType} # RNG for ionization events
    χ_v₀::ScalarField{Vector{Float64},RandomGeneratorType} # RNG for initial velocity sampling
    v₋::Vector3D{Vector{Float64},ParticleVelocityType}
    v′::Vector3D{Vector{Float64},ParticleVelocityType}
    v₊::Vector3D{Vector{Float64},ParticleVelocityType}
    s::Vector3D{Vector{Float64},ParticleVelocityType}
end

Helpers(Np::Int64) = Helpers([ft(Np) for ft in fieldtypes(Helpers)]...)

export Helpers