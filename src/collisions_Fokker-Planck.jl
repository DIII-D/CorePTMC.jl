function collision_FP!(p::ParticlePosition, v::ParticleVelocity, Z::Vector{Float64}, m::Vector{Float64}, main_ion::Element,
    nᵢ::MainIonDensity, Tᵢ::MainIonTemperature,  vᵢ::MainIonVelocity, Tₑ::ElectronTemperature, D_perp::AnomalousDiffusion, dt::Float64, i::Int64)

    mᵢ = main_ion.m # Deuterium background
    Zᵢ = main_ion.Z  # assumption of Hydrogenic plasma

    # plasma-impurity relative velocity
    wx = vᵢ.x[i] - v.x[i]
    wy = vᵢ.y[i] - v.y[i]
    wz = vᵢ.z[i] - v.z[i]
    w_mod = sqrt(wx^2 + wy^2 + wz^2)

    # I could put lnΛ in helpers
    lnΛ = 6.6 - 0.5 * log(nᵢ.value[i] / 10 ^ 20) + 1.5 * log(Tₑ.value[i]) # Te in eV and ni in m^-3; eq. 1.8 page 9 in book "The Physics of plasmas" by T.J.M Boyd and J.J. Sanderson
    # lnΛ = log(4 * π * ni * r_d^3) # classic formulation

    A = Zᵢ^2 * Z[i]^2 * ee ^ 4 / (4 * π * ϵ₀ ^ 2 * m[i] ^ 2) * lnΛ # m^6 / s^4; proportional to coulomb collision cross section. eq. 2.21 ERO-TEXTOR version 2.0 MANUAL 
    τ_ab = w_mod .^ 3 ./ ( nᵢ.value[i] * A )  # s; eq. 2.36 ERO-TEXTOR version 2.0 MANUAL
    Er = w_mod ^ 2 * mᵢ / (2 * Tᵢ.value[i] * ee) # Energy ratio between relative kinetic energy and thermal energy of background plasma

    τ_s = τ_ab / μ(Er) * mᵢ / (mᵢ + m[i]) # s; slowing down time. Eq. 2.31 ERO-TEXTOR version 2.0 MANUAL
    τ_E = Er * τ_ab / μ(Er) # s; energy exchange time. Eq. 2.31 ERO-TEXTOR version 2.0 MANUAL
    τ_d = τ_ab / ( 2 * (1 - 1 / (2 * Er) ) * μ(Er) + μ_prime(Er))  # s; deflection time. Eq. 2.31 ERO-TEXTOR version 2.0 MANUAL

    w_parallel, w_perp1, w_perp2 = orthonormal_basis_aligned_to_vector([wx, wy, wz])
    
    #@show v.x[i], v.y[i], v.z[i]
    #@show vᵢ.x[i], v.y[i], v.z[i]
    #@show w_mod
    #@show A
    #@show τ_ab
    #@show Er
    #@show τ_s, τ_E, τ_d
    
    # m/s; delta of velocity associated to Fokker-Planck operators
    Δv_FP = w_mod / τ_s * dt .* w_parallel + 
        random_direction(rand()) * sqrt(w_mod ^ 2 / τ_E * dt) .* w_parallel +
        random_direction(rand()) * sqrt(0.5 * (w_mod ^ 2 / τ_d) * dt) .* w_perp1 +
        random_direction(rand()) * sqrt(0.5 * (w_mod ^ 2 / τ_d) * dt) .* w_perp2

    # Anomalous diffusion transport 
    if Z[i] < 1.0 || D_perp.value[i] < 1e-3
        Δp_perp = [0.0, 0.0, 0.0]
    else
        Δp_perp = sqrt(D_perp.value[i] * dt) .* random_vector_in_plane(w_perp1, w_perp2)
    end

    # Update velocity
    v.x[i] = v.x[i] + Δv_FP[1]
    v.y[i] = v.y[i] + Δv_FP[2]
    v.z[i] = v.z[i] + Δv_FP[3]

    # Update position using new velocity and anomalous transport coefficient
    p.x[i] += v.x[i] * dt + Δp_perp[1]
    p.y[i] += v.y[i] * dt + Δp_perp[2]
    p.z[i] += v.z[i] * dt + Δp_perp[3]

end

## -------- HELPER FUNCTIONS -------- ##

function μ(x::Union{Float64, Vector{Float64}})   # Eq. 2.34 ERO-TEXTOR version 2.0 MANUAL
    v = 10 .^ LinRange(log10(1e-3), log10(x), 100)
    return 2 / sqrt(π) * trapz(v, exp.(-v) .* sqrt.(v) )
end

function μ_prime(x::Union{Float64, Vector{Float64}}) # Eq. 2.35 ERO-TEXTOR version 2.0 MANUAL
    return 2 / sqrt(π) .* exp.(-x) .* sqrt.(x)
end

## HELPERS FUNCTIONS

function orthonormal_basis_aligned_to_vector(w::AbstractVector{<:Real})
    # Check that the input is a 3D vector.
    if length(w) != 3
        error("Input vector must be 3-dimensional")
    end
    
    # Compute the normalized vector (w_parallel)
    norm_w = norm(w)
    if norm_w == 0
        error("Zero vector provided. Cannot compute an orthonormal basis.")
    end
    w_parallel = w / norm_w
    
    # Choose an arbitrary vector that is not parallel to w_parallel. Its cross product with vparallel will produce a vector that is perp. to w_parallel
    # We choose the coordinate axis that is least aligned with w_parallel.
    if abs(w_parallel[1]) < abs(w_parallel[2]) && abs(w_parallel[1]) < abs(w_parallel[3])
        arbitrary = [1.0, 0.0, 0.0]
    elseif abs(w_parallel[2]) < abs(w_parallel[3])
        arbitrary = [0.0, 1.0, 0.0]
    else
        arbitrary = [0.0, 0.0, 1.0]
    end
    
    # Compute w_perp1 as the cross product of w_parallel and the arbitrary vector.
    w_perp1 = cross(w_parallel, arbitrary)
    w_perp1 = w_perp1 / norm(w_perp1)  # Normalize w_perp1
    
    # Compute w_perp2 as the cross product of w_parallel and w_perp1.
    w_perp2 = cross(w_parallel, w_perp1)
    
    return w_parallel, w_perp1, w_perp2
end

# Example usage:
# w = [0.0, 1.0, 0.0]
# w_parallel, w_perp1, w_perp2 = orthonormal_basis_aligned_to_vector(w)
# println("w_parallel = ", w_parallel)
# println("w_perp1 = ", w_perp1)
# println("w_perp2 = ", w_perp2)

function random_direction(x)
    if x > 0.5
        return 1.0
    else
        return -1.0
    end
end

function random_vector_in_plane(w_perp1::AbstractVector{<:Real}, w_perp2::AbstractVector{<:Real})
    # Generate a random angle between 0 and 2π.
    theta = 2π * rand()
    # Create a random vector in the plane individuated by w_perp1 and w_perp2 using the linear combination.
    return cos(theta) * w_perp1 + sin(theta) * w_perp2
end