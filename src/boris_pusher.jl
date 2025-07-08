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
function boris_pusher!(p::ParticlePosition, v::ParticleVelocity, E::ElectricField, B::MagneticField, helpers::Helpers, dt::Float64, i::Int64)
    α = helpers.α
    v₋ = helpers.v₋
    v′ = helpers.v′
    v₊ = helpers.v₊
    s = helpers.s


    # Half acceleration due to electric field
    # α = 0.5 * (q / m) * dt 
    v₋.x[i] = v.x[i] + α[i] * E.x[i]
    v₋.y[i] = v.y[i] + α[i] * E.y[i]
    v₋.z[i] = v.z[i] + α[i] * E.z[i]

    # Rotation due to magnetic field


    v′.x[i] = v₋.x[i] + α[i] * cross_x(v₋, B, i)
    v′.y[i] = v₋.y[i] + α[i] * cross_y(v₋, B, i)
    v′.z[i] = v₋.z[i] + α[i] * cross_z(v₋, B, i)
    norm_square_B = dot(B, B, i)
    s.x[i] = 2.0 * α[i] * B.x[i] / (1.0 + α[i]^2 * norm_square_B)
    s.y[i] = 2.0 * α[i] * B.y[i] / (1.0 + α[i]^2 * norm_square_B)
    s.z[i] = 2.0 * α[i] * B.z[i] / (1.0 + α[i]^2 * norm_square_B)

    v₊.x[i] = v₋.x[i] + cross_x(v′, s, i)
    v₊.y[i] = v₋.y[i] + cross_y(v′, s, i)
    v₊.z[i] = v₋.z[i] + cross_z(v′, s, i)

    # Final update to velocity

    v.x[i] = v₊.x[i] + α[i] * E.x[i]
    v.y[i] = v₊.y[i] + α[i] * E.y[i]
    v.z[i] = v₊.z[i] + α[i] * E.z[i]

    # Update position using new velocity
    p.x[i] += v.x[i] * dt
    p.y[i] += v.y[i] * dt
    p.z[i] += v.z[i] * dt
end

update_α!(α, Z, m, dt) = @. α = 0.5 * (Z * ee / m) * dt
update_α!(α, Z, m, dt, i) = α.value[i] = 0.5 * (Z[i] * ee / m[i]) * dt