function update_particle_location!(particle::Particles, ip, dt)
    particle.p.x[ip] += particle.v.x[ip] * dt
    particle.p.y[ip] += particle.v.y[ip] * dt
    particle.p.z[ip] += particle.v.z[ip] * dt
    nothing
end

include("boris_pusher.jl")

include("collisions_Fokker-Planck.jl") # contains anomalous diffusion 

struct ForcesFlags
    boris::Bool
    FP::Bool
end

function update_position!(particle::Particles, plasma::PlasmaBackground, helpers, forces_flags::ForcesFlags, num_params::AbstractNumericalParameters, i) #helpers::PusherHelpers, dt)
    
    
    update_α!(helpers.α, particle.Z, particle.m, num_params.dt, i) # create α and store in helpers

    if forces_flags.boris && forces_flags.FP
        #println("Forces: Boris and Fokker-Planck")
        boris_pusher!(particle.p, particle.v, plasma.E, plasma.B, helpers, num_params.dt, i)
        collision_FP!(particle.p, particle.v, particle.Z, particle.m, plasma.el_bg, plasma.nᵢ, plasma.Tᵢ, plasma.vᵢ, plasma.Tₑ, plasma.D_perp, num_params.dt, i)
    elseif forces_flags.boris
        #println("Forces: Boris")
        boris_pusher!(particle.p, particle.v, plasma.E, plasma.B, helpers, num_params.dt, i) # calculate trajectories
    elseif forces_flags.FP
        #println("Forces: Fokker-Planck")
        collision_FP!(particle.p, particle.v, particle.Z, particle.m, plasma.el_bg, plasma.nᵢ, plasma.Tᵢ, plasma.vᵢ, plasma.Tₑ, plasma.D_perp, num_params.dt, i)
    else
        #@warn "No forces have been chosen by the user!"
        particle.p.x[i] += particle.v.x[i] * num_params.dt
        particle.p.y[i] += particle.v.y[i] * num_params.dt
        particle.p.z[i] += particle.v.z[i] * num_params.dt
    end

    #@show particle.p.z

    return particle
end

function update_position!(particle::Particles, plasma::PlasmaBackground, helpers, num_params::AbstractNumericalParameters, ip) #helpers::PusherHelpers, dt)

    update_α!(helpers.α, particle.Z, particle.m, num_params.dt, ip) # create α and store in helpers
    boris_velocity_pusher!(particle.v, plasma.E, plasma.B, helpers, ip) # calculate trajectories
    # if forces_flags.boris && forces_flags.FP
    #     #println("Forces: Boris and Fokker-Planck")
    #     boris_pusher!(particle.p, particle.v, plasma.E, plasma.B, helpers, num_params.dt, ip)
    #     collision_FP!(particle.p, particle.v, particle.Z, particle.m, plasma.el_bg, plasma.nᵢ, plasma.Tᵢ, plasma.vᵢ, plasma.Tₑ, plasma.D_perp, num_params.dt, ip)
    # elseif forces_flags.boris
    #     #println("Forces: Boris")
    #     boris_pusher!(particle.p, particle.v, plasma.E, plasma.B, helpers, num_params.dt, i) # calculate trajectories
    # elseif forces_flags.FP
    #     #println("Forces: Fokker-Planck")
    #     collision_FP!(particle.p, particle.v, particle.Z, particle.m, plasma.el_bg, plasma.nᵢ, plasma.Tᵢ, plasma.vᵢ, plasma.Tₑ, plasma.D_perp, num_params.dt, i)
    # else
        #@warn "No forces have been chosen by the user!"
    update_particle_location!(particle, ip, num_params.dt)
end
