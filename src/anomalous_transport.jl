struct AnomalousDiffusionTransport{D,G,T,C,E,F}
    D_perp::D
    Δt_perp::G
    dt_perp::Float64
    previous_time::T
    χ_perp::C
    e_rad::E
    e_dia::F
end


function anomalous_diffusion_pusher!(ad::AnomalousDiffusionTransport, particle_data, ip::Int64)
    particle = particle_data.particle
    particle.Z[ip] ==0 && return nothing
    particle.fly_time[ip] < ad.previous_time[ip] + ad.dt_perp && return nothing
    _anomalous_diffusion_pusher!(ad, particle, ip)
    return
end

function _anomalous_diffusion_pusher!(ad, particle, ip::Int64)

    ad.Δt_perp[ip] = particle.fly_time[ip] + particle.dt[ip] - ad.Δt_perp[ip]
    Δp = sqrt(ad.D_perp.value[i] * ad.Δt_perp[ip])
    update!(ad.χ_perp, ip)
    Δp_perp.x[ip], Δp_perp.y[ip], Δp_perp.z[ip] = random_vector_in_plane(Δp, e_rad, e_dia, χ_perp, ip)
    # Update position using new velocity and anomalous transport coefficient
    p.x[ip] += Δp_perp.x[ip]
    p.y[ip] += Δp_perp.y[ip]
    p.z[ip] += Δp_perp.z[ip]
    ad.Δt_perp[ip] = particle.fly_time[ip] + particle.dt[ip]
    return
end