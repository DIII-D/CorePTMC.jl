## -------- CONSTANTS -------- ##

const ee = 1.602176634e-19  # Elementary charge [C]
const mₑ = 9.10938356e-31  # Electron mass [kg]
const mₚ = 1.67e-27  # Proton mass [kg]
const m_amu = 1.66054e-27 # mass of 1 amu [kg]
const ϵ₀ = 8.8541878128e-12  # Vacuum permittivity [F/m]

export mₑ, mₚ, m_amu, ϵ₀, ee

## -------- BASE CUSTOM DATATYPES -------- ##

abstract type AbstractVector3D{U} end
abstract type AbstractScalarField{U} end

## Define memory structure and datatypes

# Vector3D and ScalarField objects are defined as mutable structures{T,U} of type AbstractVector3D{U} and AbstractScalarField{U} where 'U' is an addictional type used to dispatch. 
# and 'T' is the type of stored data within the structure. 

# vector structure
mutable struct Vector3D{T,U} <: AbstractVector3D{U}
    x::T
    y::T
    z::T
end


# scalar field structure
mutable struct ScalarField{T,U} <: AbstractScalarField{U} # first argument is datatype and second argument is ScalarField additional type
    value::T
end

## -------- CUSTOM IN-BUILT FUNCTIONS TO DEAL WITH CUSTOM DATATYPES -------- ##
include("utils2.jl")  # in-built functions

## -------- CUSTOM BASE DATATYPES CONSTRUCTORS -------- ## 

## constructors of memory blueprint 

# -------- Vector3D constructors -------- #

# definition of methods for AbstractVector3D and Vector3D to construct Vector3D structures

AbstractVector3D{U}(x::T, y::T, z::T) where {T,U} = Vector3D{T,U}(x, y, z)                             # constructing with user AbstractType U and three scalar inputs (x, y, z) of type T
AbstractVector3D{U}(v::Vector{T}) where {T,U} = Vector3D{T,U}(v...)                                    # constructing with vector input v = [x, y, z] 
AbstractVector3D{U}(Np::Int64) where {U} = AbstractVector3D{U}(zeros(Np), zeros(Np), zeros(Np))        # construcing Np vectors of zeros
AbstractVector3D{U}(Np::Int64, Nstore::Int64) where {U} = AbstractVector3D{U}(zeros(Np, Nstore), zeros(Np, Nstore), zeros(Np, Nstore))    # construcing Np vectors of zeros, Nstore times (option used to save the same vector at different time steps)

# default vector constructor for Float64 (Np zeros-vectors = initialization)
Vector3D{Vector{Float64},U}(Np::Int64) where {U} = AbstractVector3D{U}(Np)

# -------- ScalarField constructors -------- #

# definition of methods for AbstractScalarField and ScalarField to construct ScalarField structures

AbstractScalarField{U}(x::T) where {T,U} = ScalarField{T,U}(x)                                         # if input is a scalar of type T, creates a ScalarField structure with that value of type T
AbstractScalarField{U}(Np::Int64) where {U} = ScalarField{Vector{Float64},U}(zeros(Np))                # if input is a scalar of type Int64, creates a ScalarField structure with a 1D vector of zeros (1 scalar value for Np particles = a vector)
AbstractScalarField{U}(Np::Int64, Ny::Int64) where {U} = ScalarField{Matrix{Float64},U}(zeros(Np, Ny))  # if input is two scalars of type Int64, creates a ScalarField structure with a 1D vector of zeros Ny times
AbstractScalarField{U}(Np::Int64, Ny::Int64, Nt::Int64) where {U} = ScalarField{Array{Float64,3},U}(zeros(Np, Ny, Nt))  # if input is three scalars of type Int64, creates a ScalarField structure with a 2D matrix of zeros Nt times

# default Scalar constructor (initialization)
ScalarField{T,U}(Np::Int64) where {T,U} = AbstractScalarField{U}(Np)
ScalarField{Vector{Float64},U}(Np::Int64) where {U} = ScalarField{Vector{Float64},U}(zeros(Np))

## -------- ADDITIONAL CUSTOM DATATYPES -------- ##

# needed to call constructors from previous section. E.g. : generic_vector = Vector3D{Vector{Float64}, GenericType}(1)
abstract type GenericType end
abstract type ElectricFieldType end
abstract type MagneticFieldType end
abstract type ElectricPotentialType end
abstract type SpatialParticleCountType end
abstract type ElectronDensityType end
abstract type ElectronTemperatureType end
abstract type ParallelGradientElectronTemperatureType end
abstract type MainIonDensityType end
abstract type MainIonTemperatureType end
abstract type ParallelGradientMainIonTemperatureType end
abstract type MainIonVelocityType end
abstract type AnomalousDiffusionType end

abstract type UnitaryRadialComponentType end
abstract type UnitaryDiamagneticComponentType end
abstract type UnitaryParallelComponentType end

## Assigning types and subtypes to variables

ElectronDensity = AbstractScalarField{ElectronDensityType}
MainIonDensity = AbstractScalarField{MainIonDensityType}
ElectronTemperature = AbstractScalarField{ElectronTemperatureType}
MainIonTemperature = AbstractScalarField{MainIonTemperatureType}
ParallelGradientElectronTemperature = AbstractVector3D{ParallelGradientElectronTemperatureType}
ParallelGradientMainIonTemperature = AbstractVector3D{ParallelGradientMainIonTemperatureType}
ElectricPotential = AbstractScalarField{ElectricPotentialType}
SpatialParticleCount = AbstractScalarField{SpatialParticleCountType}
AnomalousDiffusion = AbstractScalarField{AnomalousDiffusionType}
MainIonVelocity = AbstractVector3D{MainIonVelocityType}
ElectricField = AbstractVector3D{ElectricFieldType}
MagneticField = AbstractVector3D{MagneticFieldType}
UnitaryRadialComponent = AbstractVector3D{UnitaryRadialComponentType}
UnitaryDiamagneticComponent = AbstractVector3D{UnitaryDiamagneticComponentType}
UnitaryParallelComponent = AbstractVector3D{UnitaryParallelComponentType}


abstract type AbstractNumericalParameters end
abstract type AbstractThreadedMCSimulation end
abstract type AbstractMCSimulation end
