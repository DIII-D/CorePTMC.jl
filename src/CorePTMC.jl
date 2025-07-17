module CorePTMC
using StatsBase
import StatsBase: Weights
using Random # import Random standard library
using Interpolations
using ADAS # you need ADAS package. To add it, type in REPL: ] add https://github.com/ProjectTorreyPines/ADAS.jl
# for a specific branch write in REPL: ]  https://github.com/ProjectTorreyPines/ADAS.jl#dev_branch_name
import SpecialFunctions # contains erf and erfinv

export AbstractNumericalParameters, AbstractMCSimulation, AbstractThreadedMCSimulation, AbstractSimulation,
       AbstractParticle, AbstractParticles, AbstractParticleVelocity, AbstractParticlesVelocity,
       ParticlePosition, ParticleVelocity, ParticlesPosition, ParticlesVelocity,
       ParticleProperties, ParticlesProperties, ParticleWeights, ParticlesWeights,
       SheathModel, MagneticField, PlasmaBackground,
    IonizationRates, NoIonizationRates, ParticleSurfaceDistance,
    PusherHelpers, ForcesFlags, Element, SpatialParticleCount, ScalarField, Particles, ParticlesData,
       ElectricField, MagneticField, ElectricPotential, ElectronDensity, MainIonDensity,
       ElectronTemperature, MainIonTemperature, MainIonVelocity, AnomalousDiffusion,
       IonizationRates, DeltaTimeDistribution, ExpTimeDistribution, NoIonizationRates,
    AbstractIonizationRates, IonizationRateData, Weights, AbstractParticleData, update!,
    update_ionization_state!, sputtered_velocity_distribution!, energy_thompson_distribution,
    is_redeposited,
    poloidal_cosine_distribution, azimuthal_flat_distribution, ρᵢ, nₑ_model, BoltzmannElectronDensityModel, Tₑ_model, D_perp_model
include("types.jl")
### !!!! DO NOT CHANGE ORDER OF INCLUDED FILES !!!! ###

## -------- RANDOM NUMBER GENERATOR -------- ##

include("random_number_generator.jl")

## -------- ELEMENTS DATABASE -------- ##

include("elements.jl")

## -------- PARTICLES PROPERTIES  -------- ##

include("particle_properties.jl")

## -------- PARTICLES WEIGHTS AND REDEP CRITERION -------- ##

include("particles_weights.jl")


## -------- HELPERS PARAMETERS DEFINITION  -------- ##

include("helpers.jl")

## -------- PLASMA BACKGROUND  -------- ##

include("sheath_model.jl")               # calculates sheath potential drop at particles' positions 
include("magnetic_field.jl")             # calculates magnetic field at particles' positions
include("plasma.jl")

## -------- IONIZATION RATES -------- ##

include("ionization_rates2.jl")

## -------- LPTMC MODELS FOR PARTICLES' TRANSPORT -------- ##

include("pusher.jl")                     # particles' pusher

## -------- FUNCTIONS TO SET PARTICLE SOURCE INITIAL POSITION AND VELOCITY  -------- ##

include("initial_position_distribution.jl")
include("initial_velocity_distribution.jl")

@exportAll
end
