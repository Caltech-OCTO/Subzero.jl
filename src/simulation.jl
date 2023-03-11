"""
Structs and functions to create and run a Subzero simulation
"""

@kwdef struct Constants{FT<:AbstractFloat}
    ρi::FT = 920.0              # Ice density
    ρo::FT = 1027.0             # Ocean density
    ρa::FT = 1.2                # Air density
    Cd_io::FT = 3e-3            # Ice-ocean drag coefficent
    Cd_ia::FT = 1e-3            # Ice-atmosphere drag coefficent
    Cd_ao::FT = 1.25e-3         # Atmosphere-ocean momentum drag coefficient
    f::FT = 1.4e-4              # Ocean coriolis parameter
    turnθ::FT = 15*pi/180       # Ocean turn angle
    L::FT = 2.93e5              # Latent heat of freezing [Joules/kg]
    k::FT = 2.14                # Thermal conductivity of surface ice[W/(m*K)]
    ν::FT = 0.3                 # Poisson's ratio
    μ::FT = 0.2                 # Coefficent of friction
    E::FT = 6e6                 # Young's Modulus
end


"""
    Simulation{FT<:AbstractFloat, DT<:Domain{FT}}

Simulation which holds a model and parameters needed for running the simulation.
This includes physical constants (consts), a random number generator (rng), the
number of seconds in a timestep (Δt), and the number of timesteps to run (nΔt).
We also have a flag for verbose, which will print out the number of timesteps
every 50 timesteps and a simulation name, which can be used when saving files.
The user can also define settings for each physical process.
"""
@kwdef struct Simulation{
    FT<:AbstractFloat,
    GT<:AbstractGrid,
    DT<:Domain,
    CT<:AbstractFractureCriteria,
}
    model::Model{FT, GT, DT}            # Model to simulate
    consts::Constants{FT} = Constants() # Constants used in Simulation
    rng = Xoshiro()                     # Random number generator 
    verbose::Bool = false               # String output printed during run
    name::String = "sim"                # Simulation name for printing/saving
    # Timesteps ----------------------------------------------------------------
    Δt::Int = 10                        # Simulation timestep (seconds)
    nΔt::Int = 7500                     # Total timesteps simulation runs for
    # Physical Processes -------------------------------------------------------
    coupling_settings::CouplingSettings = CouplingSettings()
    collision_settings::CollisionSettings{FT} = CollisionSettings()
    fracture_settings::FractureSettings{CT} = FractureSettings()
    simp_settings::SimplificationSettings{FT} = SimplificationSettings()
    # Output Writers -----------------------------------------------------------
    writers::OutputWriters{FT} = OutputWriters()
end

"""
timestep_sim!(sim, tstep, writers, ::Type{T} = Float64)

Run one step of the simulation and write output. 
Inputs:
    sim     <Simulation> simulation to advance
    tstep   <Int> current timestep
    writers <Vector{AbstractOutputWriter}> list of output OutputWriters
Outputs:
    None. Simulation advances by one timestep. 
"""
function timestep_sim!(sim, tstep, ::Type{T} = Float64) where T
    if !isempty(sim.model.floes)
        # Add ghost floes through periodic boundaries
        n_init_floes = length(sim.model.floes) # number of floes before ghosts
        add_ghosts!(sim.model.floes, sim.model.domain)

        # Output at given timestep
        write_data!(sim, tstep)  # Horribly type unstable
        
        # Collisions
        remove = zeros(Int, n_init_floes)
        transfer = zeros(Int, n_init_floes)
        if sim.collision_settings.collisions_on
            remove, transfer = timestep_collisions!(
                sim.model.floes,
                n_init_floes,
                sim.model.domain,
                remove,
                transfer,
                sim.consts,
                sim.Δt,
                sim.collision_settings,
                T
            )
        end

        # Remove the ghost floes - only used for collisions
        sim.model.floes = sim.model.floes[1:n_init_floes]
        empty!.(sim.model.floes.ghosts) 

        # Physical processes without ghost floes
        # Effects of ocean and atmosphere on ice and visa versa
        if sim.coupling_settings.coupling_on && mod(tstep, sim.coupling_settings.Δt) == 0
            timestep_coupling!(sim.model, sim.consts, sim.coupling_settings, T)
        end

        # Timestep ocean
        timestep_ocean!(sim.model, sim.consts, sim.Δt)
        
        # Move and update floes based on collisions and ocean/atmosphere forcing
        timestep_floe_properties!(sim.model.floes, sim.Δt)

        # Fracture floes
        if sim.fracture_settings.fractures_on && mod(tstep, sim.fracture_settings.Δt) == 0
            sim.model.max_floe_id =
                fracture_floes!(
                    sim.model.floes,
                    sim.model.max_floe_id,
                    sim.rng,
                    sim.fracture_settings,
                    sim.coupling_settings,
                    sim.simp_settings,
                    sim.consts,
                    T
                )
        end
        
        # Remove floes that were killed or are too small in this timestep
        remove_idx = findall(
            f -> !f.alive || f.area < sim.simp_settings.min_floe_area,
            sim.model.floes,
        )
        while !isempty(remove_idx)
            idx =  pop!(remove_idx)
            StructArrays.foreachfield(col -> deleteat!(col, idx), sim.model.floes)
        end
    end

    # h0 = real(sqrt.(Complex.((-2Δt * newfloe_Δt) .* hflx)))
    # mean(h0)

    #TODO: Add dissolved mass
    return 
end

"""
    run!(sim, writers, t::Type{T} = Float64)

Run given simulation and generate output for given writers.
Simulation calculations will be done with Floats of type T (Float64 of Float32).

Inputs:
    sim     <Simulation> simulation to run
    writers <Vector{:<OutputWriters}> list of output writers
    t       <Type> Float type model is running on (Float64 or Float32)
Outputs:
    None. The simulation will be run and outputs will be saved in the output
    folder. 
"""
function run!(sim, ::Type{T} = Float64) where T
    # Start simulation
    sim.verbose && println(string(sim.name, " is running!"))
    tstep = 0
    while tstep <= sim.nΔt
        if sim.verbose && mod(tstep, 50) == 0
            println(tstep, " timesteps")
        end
        # Timestep the simulation forward
        timestep_sim!(sim, tstep, T)
        tstep+=1
    end
    sim.verbose && println(string(sim.name, " done running!"))
    return
end