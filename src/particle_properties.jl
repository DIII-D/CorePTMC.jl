## -------- PARTICLE PROPERTIES TYPES -------- ##
abstract type AbstractParticleData end
abstract type ParticlePositionType end
abstract type ParticleVelocityType end
abstract type ParticleSurfaceDistanceType end
ParticlePosition = AbstractVector3D{ParticlePositionType}
ParticleVelocity = AbstractVector3D{ParticleVelocityType}
ParticleSurfaceDistance = AbstractScalarField{ParticleSurfaceDistanceType}

## -------- PARTICLE STRUCTURE AND CONSTRUCTOR  -------- ##

# contains information about transported particles  

struct Particles{P1,P2,V,D, W ,T,U}
    p::P1                # Position (x, y, z)
    v::V                 # Velocity (vx, vy, vz)
    δsurf::D             # Distance from closest surface [m]
    w::W                 # Weights
    Z::T                 # Charge state
    m::T                 # Mass [kg]
    Zmax::T              # Atomic number
    t_fly::T             # Flying time [s]
    p_iz::P2             # Position of ionization 
    t_iz::T              # Time since last ionization
    redeposited::U       # Is particle redeposited
    id::U                # Particle ID
end

# initializing Particles structure using Element properties for a number Np of simulated particles (Element structure defined in elements.jl)
Particles(el::Element, Np::Int64; kw...) = Particles(el.Z::Float64, el.m::Float64, el.Zmax::Float64, Np::Int64; kw...)


function Particles(Z::Float64, m::Float64, Zmax::Float64, Np::Int64; ids::Union{Vector{Int},Missing}=missing)

    @assert Np > 0 "Number of particles (Np) must be greater than 0"
    if ids isa Missing
        ids = collect(1:Np)
    else
        @assert length(ids) == Np "particles id provided do not match the number of particles: np = $Np;  length(ids) = $(length(ids))"
    end
    # Debugging: Print out types
    println("Creating Particles with Z = $Z (", typeof(Z), "), m = $m (", typeof(m), "), Zmax = $Zmax (", typeof(Zmax), "), Np = $Np (", typeof(Np), ") ")

    # constructing particle structure
    return Particles(
        ParticlePosition(Np),
        ParticleVelocity(Np),
        ParticleSurfaceDistance(Np),
        weights(fill(1.0, Np)),
        Vector{Float64}(fill(Z, Np)),
        Vector{Float64}(fill(m, Np)),
        Vector{Float64}(fill(Zmax, Np)),
        Vector{Float64}(fill(0.0, Np)),
        ParticlePosition(Np),
        Vector{Float64}(fill(0.0, Np)),
        Vector{Int}(fill(0, Np)),
        ids #Vector{Int}(collect(1:Np))
    )
end