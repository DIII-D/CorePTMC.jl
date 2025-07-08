## -------- EXTRA FUNCTIONS TO PROCESS I/O MC SIMULATION  -------- ##

# set initial position of particle source !!!

# in mc_simulation the code stores initialized variables and parameters in different structures that are than used to start the simulation
# with function run! in run.jl
# you can modify the values in these structure before they are passed in run! to set custom initial velocity and position

# For instance: structures are stored in variable sim of type MCSimulation (sim = mc_simulation(...)) with all structure that than you can modify

set_position!(sim::AbstractMCSimulation, args...; kw...) = set_position!(sim.data.particles.p, args...; kw...)
set_position!(sims::AbstractThreadedMCSimulation, args...; kw...) = [set_position!(sim.data.particles.p, args...; kw...) for sim in sims.threaded_sims]

function set_position!(p::ParticlePosition, p0::Vector{Float64}; kw...)
    p.x .= p0[1]
    p.y .= p0[2]
    p.z .= p0[3]
end

set_position!(p::ParticlePosition, p0::Float64; kw...) = set_position!(p, [p0, p0, p0]; kw...)