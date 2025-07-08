set_weight!(p::Particles; kw...) = error() #TODO

# redeposition flag function
is_redeposited(ps, i) = ps.redeposited[i] == 1

# if particle's z-position is below z_surface then, particles.redeposited = 1 (particle is redeposited) otherwise: nothing 
update_redeposited!(particles, enforce_redeposition, i, ; z_surface=0.0) = enforce_redeposition && particles.p.z[i] <= z_surface ? particles.redeposited[i] = 1 : nothing
