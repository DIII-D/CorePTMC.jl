#=
Author: Jerome Guterl (guterlj@fusion.gat.com), Luca Cappelli (cappellil@fusion.gat.com)
Company: General Atomics
ionization_rates.jl (c) 2024=#



abstract type AbstractIonizationRates{M} end

## -------- IONIZATION RATES STRUCTURE AND CONSTRUCTOR  -------- ##
# Structure for Ionizations
struct IonizationRates{D,T,U,M} <: AbstractIonizationRates{M}
    data::D  # source data of SCD [m^3 s^-1] (Effective ionization coefficients)
    ν_iz::T  # ionization frequency [s^-1]
    f_iz::U  # ionization factor [#] (used to make parametric studies)
end

struct NoIonizationRates <: AbstractIonizationRates{Missing} end  # in case ionizations are not simulated

# choose here your source data for ionization coefficients (creation of iz_rates)
adas_ionization_rates(el::Element; kw...) = ADAS.get_ionization_rate(el.symbol; kw...) # renaming function and extrapolate data in a structure

import ADAS:get_ionization_rate

adas_combined_ionization_rates(el::Element, extra_data; kw...) = ADAS.get_ionization_rate(el.symbol, extra_data; kw...) # renaming function and extrapolate data in a structure


function ADAS.get_ionization_rate(imp::Union{String,Symbol}, extra_data; kw...)
    scd = retrieve_ADAS_data(imp; type="scd", kw...)

    ndens = length(scd.data.axis.ne)
    ntemp = length(scd.data.axis.Te)
    nZ = length(scd.data.rates)
    rate = zeros(nZ, ndens, ntemp)
    #@show scd
    for Z in 1:nZ
        rate[Z, :, :] .= scd.data.rates[Z].values[:, :]
    end
    println("rate = $(rate[1, 10, 10])")

    dens, temp = ADAS.meshgrid(scd.data.axis.ne, scd.data.axis.Te)

    for Z in extra_data.Z
        println("Z extra = $Z")
        rate[Z + 1, :, :] .= extra_data.rate.(Z .+ 0 .* dens , dens, temp)
    end
    
    println("rate = $(rate[1, 10, 10])")

    Te = scd.data.axis.Te
    ne = scd.data.axis.ne
    Z = collect(0:nZ-1)
    rate_ = Interpolations.linear_interpolation((float.(Z), ne, Te), rate; extrapolation_bc=Interpolations.Flat())

    return ADAS.IonizationRate(scd, Te, ne, Z, rate_, imp)
end


no_ionization_rates() = NoIonizationRates()

# ionization time distributions
abstract type TimeDistribution end
struct ExpTimeDistribution <: TimeDistribution end    # exponential time distribution (probability of ionization increases exponentially with time, use of random generator to get ioniz. prob at each time step)
struct DeltaTimeDistribution <: TimeDistribution end  # delta time distribution (if t > t_iz particle undergoes ionization otherwise it doesn't)

# choose here your ionization time distribution
function IonizationRates(Np::Int64, iz_rates; time_distr = ExpTimeDistribution, f_iz::Float64=1.0, kw...) # constructor function for IonizationRates struct
    AbstractIonizationRates{time_distr}(iz_rates, zeros(Np), f_iz)
end

# No Ionizations case
function IonizationRates(Np::Int64, iz_rates::NoIonizationRates; time_distr = ExpTimeDistribution, f_iz::Float64=1.0, kw...) # constructor function for IonizationRates struct
    NoIonizationRates()
end

# rule to determine ionization envent
ionization_event(χ::RandomGenerator, iz::AbstractIonizationRates{ExpTimeDistribution}, dt::Float64, i) = heaviside((1.0 - exp(-iz.f_iz * iz.ν_iz[i] * dt)) - χ[i])
ionization_event(iz::AbstractIonizationRates{DeltaTimeDistribution}, t::Float64, i) = heaviside(t - 1 / iz.ν_iz[i]) # t is time from previous ionization

# assign method to AbstractIonizationRates type --> method is for collecting data into IonizationRates struct
AbstractIonizationRates{M}(data::D, ν_iz::T, f_iz::U) where {M,D,T,U} = IonizationRates{D,T,U,M}(data, ν_iz, f_iz)

# create function for no_ionization_rate case, that is for an input that is of type NoIonizationRates. In this case: do nothing let ps.Z = ps.Z at t=0
update_ionization_state!(χ::RandomGenerator, iz::NoIonizationRates, ps::Particles, bg::PlasmaBackground, dt, i) = nothing


function update_ionization_state!(χ::RandomGenerator, iz::IonizationRates, ps::Particles, bg::PlasmaBackground, dt::Float64, i)

    # @show 1.0 - exp(-iz.f_iz * iz.ν_iz[i] * ps.t_iz[i]), χ[i]

    if ps.Z[i] < ps.Zmax[i]
        update_iz_rates!(iz, ps.Z, bg.nₑ, bg.Tₑ, i)
        if ionization_event(χ, iz, dt, i) == 1.0
            ps.Z[i] += 1
            ps.t_iz[i] = 0.0
            ps.p_iz.x[i] = ps.p.x[i]
            ps.p_iz.y[i] = ps.p.y[i]
            ps.p_iz.z[i] = ps.p.z[i]
        end
    end
end

function update_ionization_state!(iz::IonizationRates{<:Any, <:Any, <:Any, <:DeltaTimeDistribution}, ps::Particles, bg::PlasmaBackground, t::Float64, i)

    if ps.Z[i] < ps.Zmax[i]
        update_iz_rates!(iz, ps.Z, bg.nₑ, bg.Tₑ, i)
        if ionization_event(iz, t, i) == 1.0
            ps.Z[i] += 1
            ps.t_iz[i] = 0.0
            ps.p_iz.x[i] = ps.p.x[i]
            ps.p_iz.y[i] = ps.p.y[i]
            ps.p_iz.z[i] = ps.p.z[i]
        end
    end
end

#update_ionization_state!(χ, iz::NoIonizationRates, ps, bg, i) = nothing

# iz.data.rate(Z[i], ne.value[i], Te.value[i]) contains SCD values in [m^3 s^-1] from ADAS package stored in IonizationRates structure
#
# SCD * nₑ is equal to the ionization rate ν_iz  

update_iz_rates!(iz::IonizationRates, Z, ne, Te, i) = iz.ν_iz[i] = iz.data.rate(Z[i], ne.value[i], Te.value[i]) * ne.value[i]