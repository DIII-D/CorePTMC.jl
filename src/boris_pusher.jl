#vectorized version


function boris_pusher!(p::ParticlePosition, v::ParticleVelocity, E::ElectricField, B::MagneticField, helpers::Helpers, dt::Float64)
    α = helpers.α
    v₋ = helpers.v₋
    v′ = helpers.v′
    v₊ = helpers.v₊
    s = helpers.s

    # Half acceleration due to electric field (charge to mass ratio)
    # α = 0.5 * (q / m) * dt 
    @. v₋.x = v.x + α * E.x
    @. v₋.y = v.y + α * E.y
    @. v₋.z = v.z + α * E.z

    # Rotation due to magnetic field

    @. v′.x = v₋.x + α * cross_x(v₋, B)
    @. v′.y = v₋.y + α * cross_y(v₋, B)
    @. v′.z = v₋.z + α * cross_z(v₋, B)

    norm_square_B = dot(B, B)

    @. s.x = 2.0 * α * B.x / (1.0 + α^2 * norm_square_B)
    @. s.y = 2.0 * α * B.y / (1.0 + α^2 * norm_square_B)
    @. s.z = 2.0 * α * B.z / (1.0 + α^2 * norm_square_B)

    @. v₊.x = v₋.x + cross_x(v′, s)
    @. v₊.y = v₋.y + cross_y(v′, s)
    @. v₊.z = v₋.z + cross_z(v′, s)

    # Final update to velocity

    @. v.x = v₊.x + α * E.x
    @. v.y = v₊.y + α * E.y
    @. v.z = v₊.z + α * E.z

    # Update position using new velocity
    @. p.x += v.x * dt
    @. p.y += v.y * dt
    @. p.z += v.z * dt
end

#kernel version
function boris_pusher!(p::ParticlePosition, v::ParticleVelocity, E::ElectricField, B::MagneticField, helpers::Helpers, dt::Float64, ip::Int64)
    α = helpers.α
    v₋ = helpers.v₋
    v′ = helpers.v′
    v₊ = helpers.v₊
    s = helpers.s


    # Half acceleration due to electric field
    # α = 0.5 * (q / m) * dt
    v₋.x[ip] = v.x[ip] + α[ip] * E.x[ip]
    v₋.y[ip] = v.y[ip] + α[ip] * E.y[ip]
    v₋.z[ip] = v.z[ip] + α[ip] * E.z[ip]

    # Rotation due to magnetic field


    v′.x[ip] = v₋.x[ip] + α[ip] * cross_x(v₋, B, ip)
    v′.y[ip] = v₋.y[ip] + α[ip] * cross_y(v₋, B, ip)
    v′.z[ip] = v₋.z[ip] + α[ip] * cross_z(v₋, B, ip)
    norm_square_B = dot(B, B, ip)
    s.x[ip] = 2.0 * α[ip] * B.x[ip] / (1.0 + α[ip]^2 * norm_square_B)
    s.y[ip] = 2.0 * α[ip] * B.y[ip] / (1.0 + α[ip]^2 * norm_square_B)
    s.z[ip] = 2.0 * α[ip] * B.z[ip] / (1.0 + α[ip]^2 * norm_square_B)

    v₊.x[ip] = v₋.x[ip] + cross_x(v′, s, ip)
    v₊.y[ip] = v₋.y[ip] + cross_y(v′, s, ip)
    v₊.z[ip] = v₋.z[ip] + cross_z(v′, s, ip)

    # Final update to velocity

    v.x[ip] = v₊.x[ip] + α[ip] * E.x[ip]
    v.y[ip] = v₊.y[ip] + α[ip] * E.y[ip]
    v.z[ip] = v₊.z[ip] + α[ip] * E.z[ip]

    # Update position using new velocity
    p.x[ip] += v.x[ip] * dt
    p.y[ip] += v.y[ip] * dt
    p.z[ip] += v.z[ip] * dt
end

#kernel version
function boris_velocity_pusher!(v::ParticleVelocity, E::ElectricField, B::MagneticField, helpers::Helpers, ip::Int64)
    α = helpers.α
    v₋ = helpers.v₋
    v′ = helpers.v′
    v₊ = helpers.v₊
    s = helpers.s


    # Half acceleration due to electric field
    # α = 0.5 * (q / m) * dt
    v₋.x[ip] = v.x[ip] + α[ip] * E.x[ip]
    v₋.y[ip] = v.y[ip] + α[ip] * E.y[ip]
    v₋.z[ip] = v.z[ip] + α[ip] * E.z[ip]

    # Rotation due to magnetic field


    v′.x[ip] = v₋.x[ip] + α[ip] * cross_x(v₋, B, ip)
    v′.y[ip] = v₋.y[ip] + α[ip] * cross_y(v₋, B, ip)
    v′.z[ip] = v₋.z[ip] + α[ip] * cross_z(v₋, B, ip) 
    γ = 1.0 + α[ip]^2 * dot(B, B, ip)
    s.x[ip] = 2.0 * α[ip] * B.x[ip] / γ
    s.y[ip] = 2.0 * α[ip] * B.y[ip] / γ
    s.z[ip] = 2.0 * α[ip] * B.z[ip] / γ

    v₊.x[ip] = v₋.x[ip] + cross_x(v′, s, ip)
    v₊.y[ip] = v₋.y[ip] + cross_y(v′, s, ip)
    v₊.z[ip] = v₋.z[ip] + cross_z(v′, s, ip)

    # Final update to velocity

    v.x[ip] = v₊.x[ip] + α[ip] * E.x[ip]
    v.y[ip] = v₊.y[ip] + α[ip] * E.y[ip]
    v.z[ip] = v₊.z[ip] + α[ip] * E.z[ip]
end

update_α!(α, Z, m, dt) = @. α = 0.5 * (Z * ee / m) * dt
update_α!(α, Z, m, dt, i) = α.value[i] = 0.5 * (Z[i] * ee / m[i]) * dt