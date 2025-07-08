#=
Authors: Jerome Guterl (guterlj@fusion.gat.com), Luca Cappelli (cappellil@fusion.gat.com)
 Company: General Atomics
 velocity_distribution.jl (c) 2024=#

# modify initial velocity data stored in structure
set_velocity!(sim::AbstractMCSimulation, args...; kw...) = set_velocity!(sim.data.particles.v, args...; kw...)
set_velocity!(sims::AbstractThreadedMCSimulation, args...; kw...) = [set_velocity!(sim.data.particles.v, args...; kw...) for sim in sims.threaded_sims]
## delta distributions: all particles simulated with the same initial velocity
function set_velocity!(v::ParticleVelocity, v0::Vector{Float64}; kw...)
    v.x .= v0[1]
    v.y .= v0[2]
    v.z .= v0[3]
end
# methods of function set_velocity!
set_velocity!(v::ParticleVelocity, v0::Float64; kw...) = set_velocity!(v, [v0, v0, v0]; kw...)

## flat distribution: particles with velocity ∈ [a, b]

# modify initial velocity data stored in structure
velocity_flat_distribution!(sim::AbstractMCSimulation, args...; kw...) = velocity_flat_distribution!(sim.data.particles.v, args...; kw...)

function velocity_flat_distribution!(v::ParticleVelocity, a::Vector{Float64}, b::Vector{Float64})
    Np = length(v.x)
    χx = random_vector(Np)
    χy = random_vector(Np)
    χz = random_vector(Np)
    Δ = b - a
    @. v.x = χx.value * Δ[1] + a[1]
    @. v.y = χy.value * Δ[2] + a[2]
    @. v.z = χz.value * Δ[3] + a[3]
end


## Maxwellian distribution sampling

velocity_maxwellian_distribution!(sim::AbstractMCSimulation, args...; kw...) = velocity_maxwellian_distribution!(sim.data.particles.v, args...; kw...)
velocity_maxwellian_distribution!(sims::AbstractThreadedMCSimulation, args...; kw...) = [velocity_maxwellian_distribution!(sim.data.particles.v, args...; kw...) for sim in sims.threaded_sims]

function velocity_maxwellian_distribution!(v::ParticleVelocity, c::Vector{Float64}, T::Float64, el::Element)

    # c = average velocity vector
    # T = temperature

    Np = length(v.x)
    β = el.m / (2 * T * ee)
    χx = random_vector(Np)
    χy = random_vector(Np)
    χz = random_vector(Np)
    v.x = @. c[1] + erfinv(2 * χx.value - 1) / sqrt(β)
    @. v.y = c[2] + erfinv(2 * χy.value - 1) / sqrt(β)
    @. v.z = c[3] + erfinv(2 * χz.value - 1) / sqrt(β)
end

## 1D Maxwellian distribution

function maxwellian_1D_distribution(v::Vector{Float64}, c::Float64, T::Float64, el::Element)
    β = el.m / (2 * T * ee)
    f = @. sqrt(β / pi) * exp(-β * (v - c)^2)
    return f
end


## Any distribution sampling
# calculate Cumulative Density Function numerically for functions hard or impossible to integrate analytically
#
# hypothesis: independent variables, CDF can be calculated for each dimension


function distribution_sampling(x::Vector{Float64}, PDF::Vector{Float64}, Np::Union{Int64,Float64}; seed=123)

    CDF = cumtrapz(x, PDF)
    itp = linear_interpolation(Interpolations.deduplicate_knots!(CDF), x)

    χ = random_vector(Np)

    x = itp(χ.value)
end

## Energy distribution sampling

# particles are emitted according to spherical coordinate system where:
# E: energy
# θ: poloidal angle    ∈ [0, π/2] (does not include angles until -π/2 because this function describes sputtered particles)
# ϕ: azimuthal angle   ∈ [0, 2π]
# dvₓdvydvz = sin(θ)dθdϕ 
# PDF(E, θ, ϕ) = PDF(E) * PDF(θ) * PDF(ϕ) * sin(θ)dEdθdϕ

# E₀: Average incident energy

# Thompson

logrange(a, b, c) = 10 .^ LinRange(log10(a), log10(b), c)
function energy_thompson_distribution(E₀::Float64, el_projectile::Element, el_target::Element, Np::Int64)
    
    Ecutoff = calculate_Ecutoff(el_target, el_projectile, E₀)
    return energy_thompson_distribution(Ecutoff, el_target, Np)
end
    
function energy_thompson_distribution(Ecutoff::Float64, el_target::Element, Np::Int64)
    Es = el_target.Es
    E = collect(logrange(0.01, Ecutoff * 2, 50))
    fE = @. E / (E + Es)^3 * (1 - sqrt((E + Es) / (Ecutoff + Es))) * (E <= Ecutoff)

    if isnan(Es)
        error("Es for $(String(el_target.symbol)) not present in database, update it in src/elements.jl")
    end


    PDF = fE ./ trapz(E, fE)

    E_sampled = distribution_sampling(E, PDF, Np)

    return E_sampled, E, PDF
end

function weighted_energy_thompson_distribution(E₀::Vector{Float64}, el_projectile::Vector{Element}, w::Vector{Float64}, el_target::Element, Np::Int64; NE::Int64 = 50)
    
    Ecutoff = [calculate_Ecutoff(el_target, el_projectileᵢ, E₀ᵢ) for (E₀ᵢ, el_projectileᵢ) in zip(E₀, el_projectile)]
    return weighted_energy_thompson_distribution(Ecutoff, w, el_target, Np; NE = NE)
end

function weighted_energy_thompson_distribution(Ecutoff::Vector{Float64}, w::Vector{Float64}, el_target::Element, Np::Int64; NE::Int64 = 50)
    ϵ = 1e-6
    @assert (1 - ϵ < sum(w)) & (sum(w) < 1 + ϵ) "Weights sum must be equal to 1"
    @assert length(Ecutoff) == length(w) "length of Ecutoff must be equal to the number of weights"

    Es = el_target.Es
    E = collect(logrange(0.01, maximum(Ecutoff) * 2, NE))
    fE = zeros(NE)
    
    for (wi, Ecutoffi) in zip(w, Ecutoff) 
        fE .+= @. wi * (E / (E + Es)^3 * (1 - sqrt((E + Es) / (Ecutoffi + Es))) * (E <= Ecutoffi))
    end
    PDF = fE ./ trapz(E, fE)

    E_sampled = distribution_sampling(E, PDF, Np)

    return E_sampled, E, PDF
end

function calculate_Ecutoff(el_target, el_projectile, E₀)
    
    mt = el_target.m
    mp = el_projectile.m

    g = 4 * mp * mt / (mp + mt)^2

    return E₀ .* g * ((mp + 2 * mt) / (2 * mp + 2 * mt))^6
end

## Poloidal distribution sampling

# cosine
function poloidal_cosine_distribution(Np::Int64)

    χ = random_vector(Np)

    # In 3D poloidal cosine distribution = cos(θ)sin(θ) if normalized to make a PDF = sin(2 .* θ)
    # Inversion method: x = 0.5 .* acos(1 - 2 .* χ) 

    θ_sampled = acos.(1 .- 2 .* χ.value) .* 0.5
    return θ_sampled
end

#cosine 2D
function poloidal_cosine_distribution2D(Np)
    χ = random_vector(Np)

    # 2D poloidal cosine distribution = cos(θ) if normalized to make a PDF = 0.5 * cos(θ); CDF = 0.5*(sin(θ) + 1)
    # Inversion method: θ = asin(2 .* χ - 1) 

    θ_sampled = asin.(2 .* χ.value .- 1)
    return θ_sampled
end

## Azimuthal distribution sampling

# flat
function azimuthal_flat_distribution(Np::Int64)
    χ = random_vector(Np)
    Δ = 2 * pi
    ϕ_sampled = χ.value .* Δ
    return ϕ_sampled
end

# two values function
function azimuthal_two_values_distribution(Np::Int64, ϕ₀::Float64, ϕ₁::Float64)
    χ = random_vector(Np)
    ϕ_sampled = zeros(Np)

    ϕ_sampled[χ.value .< 0.5] .= ϕ₀
    ϕ_sampled[χ.value .>= 0.5] .= ϕ₁

    return ϕ_sampled
end


## converting from (E, θ, ϕ) system of coordinates to (vₓ vy, vz) 

function convert_2_velocity!(v::ParticleVelocity, E_sampled::Vector{Float64}, θ_sampled::Vector{Float64}, ϕ_sampled::Vector{Float64}, el::Element)
    v_mod = sqrt.(2 * ee .* E_sampled ./ el.m)
    @. v.x = v_mod * sin(θ_sampled) * cos(ϕ_sampled)
    @. v.y = v_mod * sin(θ_sampled) * sin(ϕ_sampled)
    @. v.z = v_mod * cos(θ_sampled)
end

## converting from (v, θ) system of coordinates to (vx = 0, vy = v*sin(θ), vz = v*cos(θ)) 

function convert_2_velocity!(v::ParticleVelocity, E_sampled::Vector{Float64}, θ_sampled::Vector{Float64}, el::Element)
    v_mod = sqrt.(2 * ee .* E_sampled ./ el.m)

    println("2D cosine poloidal sampling.")

    @. v.x = 0.0
    @. v.y = v_mod * sin(θ_sampled)
    @. v.z = v_mod * cos(θ_sampled)

    println("vy = $(v.y)")
end

## calculate a sputtered velocity distribution

sputtered_velocity_distribution!(sim::AbstractMCSimulation, args...; kw...) = sputtered_velocity_distribution!(sim.data.particles.v, args...; kw...)
sputtered_velocity_distribution!(sims::AbstractThreadedMCSimulation, args...; kw...) = [sputtered_velocity_distribution!(sim.data.particles.v, args...; kw...) for sim in sims.threaded_sims]

function sputtered_velocity_distribution!(v::ParticleVelocity, el_projectile::Element, el_target::Element, E₀::Float64,
    fE::Function,
    fθ::Function,
    fϕ::Function)
    Np = length(v.x)
    E_sampled, _, _ = fE(E₀, el_projectile, el_target, Np)
    θ_sampled = fθ(Np)
    ϕ_sampled = fϕ(Np)

    convert_2_velocity!(v, E_sampled, θ_sampled, ϕ_sampled, el_target)

end

function sputtered_velocity_distribution!(v::ParticleVelocity, el_target, Ecutoff,
    fE,
    fθ,
    fϕ)
    Np = length(v.x)
    E_sampled, _, _ = fE(Ecutoff, el_target, Np)
    θ_sampled = fθ(Np)
    ϕ_sampled = fϕ(Np)

    convert_2_velocity!(v, E_sampled, θ_sampled, ϕ_sampled, el_target)

end

# weighted energy distributions

function sputtered_velocity_distribution!(v::ParticleVelocity, el_projectile::Vector{Element}, el_target::Element, E₀::Vector{Float64}, w::Vector{Float64},
    fE,
    fθ,
    fϕ)
    Np = length(v.x)
    E_sampled, _, _ = fE(E₀, el_projectile, w, el_target, Np)
    θ_sampled = fθ(Np)
    ϕ_sampled = fϕ(Np)

    convert_2_velocity!(v, E_sampled, θ_sampled, ϕ_sampled, el_target)

end

function sputtered_velocity_distribution!(v::ParticleVelocity, el_target::Element, Ecutoff::Vector{Float64}, w::Vector{Float64},
    fE,
    fθ,
    fϕ)
    Np = length(v.x)
    E_sampled, _, _ = fE(Ecutoff, w, el_target, Np)
    θ_sampled = fθ(Np)
    ϕ_sampled = fϕ(Np)

    convert_2_velocity!(v, E_sampled, θ_sampled, ϕ_sampled, el_target)

end
## DISTRIBUTION TO BENCHMARK WITH FUSSMAN MODEL Fussman G et al 1995 PPCF Research 1994 Proc. 15th Int. Conf. (Seville, 1994) ed Vienna IAEA vol 2 pp 143–8 IAEA, Vienna

function sputtered_velocity_distribution!(v::ParticleVelocity, el_target, Ecutoff,
    fE,
    fθ)
    Np = length(v.x)
    E_sampled, _, _ = fE(Ecutoff, el_target, Np)
    θ_sampled = fθ(Np)

    println("2D velocity distribution.")
    convert_2_velocity!(v, E_sampled, θ_sampled, el_target)

end

#=
function sputtered_velocity_distribution!(v::ParticleVelocity, el_projectile, el_target, E₀,
    fE,
    fθ)
    Np = length(v.x)
    E_sampled, _, _ = fE(E₀, el_projectile, el_target, Np)
    θ_sampled = fθ(Np)

    convert_2_velocity!(v, E_sampled, θ_sampled, el_target)

end
=#