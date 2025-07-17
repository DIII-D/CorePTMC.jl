abstract type AbstractPlasmaModel end

struct FieldAlignedUnitVector{T,U,V}
    e_rad::T
    e_dia::U
    e_para::V
end

FieldAlignedUnitVector(Np::Int64) = FieldAlignedUnitVector(
    UnitaryRadialComponent(Np),
    UnitaryDiamagneticComponent(Np),
    UnitaryParallelComponent(Np)
)


vcat_preserve_type(vec::Vector{<:AbstractVector3D}) = vcat_vector3D(vec)
function vcat_preserve_type(vec::Vector{<:FieldAlignedUnitVector})
    return FieldAlignedUnitVector([vcat_preserve_type([getproperty(v,fn) for v in vec]) for fn in fieldnames(FieldAlignedUnitVector)]...)
end
struct PlasmaBackground{Bg,Ne,Ni,Te,Ti,GTe,GTi,Vi,U,V,W,DP,F}  # background data
    el_bg::Bg
    nₑ::Ne
    nᵢ::Ni    
    Tₑ::Te
    Tᵢ::Ti
    ∇Tₑ::GTe
    ∇Tᵢ::GTi
    vᵢ::Vi  # plasma main ion average velocity
    E::U
    B::V
    ϕ::W
    D_perp::DP  # anomalous transport
    e::F
end

# initialization of plasma background

PlasmaBackground(Np::Int64, main_ion::Element) = PlasmaBackground(
    main_ion,
    ElectronDensity(Np),        # constructors defined in 'LPTMC.jl'
    MainIonDensity(Np),
    ElectronTemperature(Np),
    MainIonTemperature(Np), 
    ParallelGradientElectronTemperature(Np),
    ParallelGradientMainIonTemperature(Np),
    MainIonVelocity(Np),
    ElectricField(Np), 
    MagneticField(Np), 
    ElectricPotential(Np),
    AnomalousDiffusion(Np),
    FieldAlignedUnitVector(Np)
)

PlasmaBackground(Np::Int64) = PlasmaBackground(
    missing,
    ElectronDensity(Np),        # constructors defined in 'LPTMC.jl'
    MainIonDensity(Np),
    ElectronTemperature(Np),
    MainIonTemperature(Np),
    ParallelGradientElectronTemperature(Np),
    ParallelGradientMainIonTemperature(Np),
    MainIonVelocity(Np),
    ElectricField(Np),
    MagneticField(Np),
    ElectricPotential(Np),
    AnomalousDiffusion(Np),
    FieldAlignedUnitVector(Np)
)

PlasmaBackground(Nc::Int64, Nv::Int64) = PlasmaBackground(
    missing,
    ElectronDensity(Nc, Nv),        # constructors defined in 'LPTMC.jl'
    MainIonDensity(Nc, Nv),
    ElectronTemperature(Nc, Nv),
    MainIonTemperature(Nc, Nv),
    ParallelGradientElectronTemperature(Nc, Nv),
    ParallelGradientMainIonTemperature(Nc, Nv),
    MainIonVelocity(Nc, Nv),
    ElectricField(Nc, Nv),
    MagneticField(Nc, Nv),
    ElectricPotential(Nc, Nv),
    AnomalousDiffusion(Nc, Nv),
    FieldAlignedUnitVector(Nc, Nv)
)



## -------- PLASMA BACKGROUND ELECTRON DENSITY MODELS -------- ##

include("electron_density.jl")

#TODO --- Luca: we could add a kinetic tail to the electron distribution and have f_e instead of n_e but this is for mid-term/long-term future development ---  

## -------- PLASMA BACKGROUND ELECTRON TEMPERATURE MODELS -------- ##

include("electron_temperature.jl")

## -------- PLASMA BACKGROUND ION DENSITY MODELS -------- ##

include("mainion_density.jl")

## -------- PLASMA BACKGROUND ION TEMPERATURE MODELS -------- ##

include("mainion_temperature.jl")

## -------- PLASMA BACKGROUND ION AVERAGE VELOCITY -------- ##

include("mainion_velocity.jl")

## -------- PLASMA BACKGROUND ION AVERAGE VELOCITY -------- ##

include("dperp.jl")


## -------- BACKGROUND MODELS STRUCTURE (PLASMA + MAGNETIC FIELD) -------- ##

struct BackgroundModels{BB,S,PED,PET,PID,PIT,PIV,DP} # background models
    B::BB                              # magnetic field model
    sheath::S                         # sheath model
    plasma_elec_dens::PED             # plasma electron density background model (e.g. Boltzmann factor)
    plasma_elec_temp::PET             # plasma electron temp background model (e.g. constant)
    plasma_mainion_dens::PID          # plasma main ion density background model (e.g. Boltzmann factor)
    plasma_mainion_temp::PIT          # plasma main ion temp background model (e.g. constant)
    plasma_mainion_vel::PIV           # plasma main ion velocity background model (e.g. constant)
    plasma_dperp::DP                  # plasma anomalous transport background model (e.g. constant)
end


# calculates plasma quantities at particle's position
update_plasma!(plasma, particles, models::BackgroundModels, ip::Int64) = update_plasma!(
    plasma, particles, 
    models.plasma_elec_dens, 
    models.plasma_elec_temp,
    models.plasma_mainion_dens,
    models.plasma_mainion_temp,
    models.plasma_mainion_vel,
    models.plasma_dperp,
    models.sheath, 
    models.B, 
    ip
)

function update_plasma!(plasma, particles, elec_dens_model::AbstractPlasmaModel, elec_temp_model::AbstractPlasmaModel,
    mainion_dens_model::AbstractPlasmaModel, mainion_temp_model::AbstractPlasmaModel, mainion_vel_model::AbstractPlasmaModel,
    dperp_model::AbstractPlasmaModel, sheath_model::AbstractSheathModel, B_model::Union{AbstractMagneticFieldModel,MagneticField}, ip::Int64)

    update_ϕ!(plasma.ϕ, particles.p, sheath_model, ip)                       # function update_ϕ! in 'sheath_model.jl'
    update_E!(plasma.E, particles.p, sheath_model, ip)                       # function update_E! in 'sheath_model.jl'
    update_B!(plasma.B, particles.p, B_model, ip)                            # function update_B! in 'magnetic_field.jl'
    update_nₑ!(plasma.nₑ, particles.p, plasma.ϕ, elec_dens_model, ip)        # function update_nₑ! in 'electron_density.jl'
    update_Tₑ!(plasma.Tₑ, particles.p, elec_temp_model, ip)                  # function update_Tₑ! in 'electron_temperature.jl'
    update_nᵢ!(plasma.nᵢ, particles.p, plasma.ϕ, mainion_dens_model, ip)     # function update_nᵢ! in 'mainion_density.jl'
    update_Tᵢ!(plasma.Tᵢ, particles.p, mainion_temp_model, ip)               # function update_Tᵢ! in 'mainion_temperature.jl'
    update_vᵢ!(plasma.vᵢ, particles.p, mainion_vel_model, ip)                # function update_vᵢ! in 'mainion_velocity.jl'
    update_dperp!(plasma.D_perp, particles.p, dperp_model, ip)               # function update_dperp! in 'dperp.jl'
end