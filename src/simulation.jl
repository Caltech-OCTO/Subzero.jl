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
    f::FT = 1.4e-4              # Ocean coriolis frequency
    turnθ::FT = 15π/180         # Ocean turn angle
    L::FT = 2.93e5              # Latent heat of freezing [Joules/kg]
    k::FT = 2.14                # Thermal conductivity of surface ice[W/(m*K)]
    ν::FT = 0.3                 # Poisson's ratio
    μ::FT = 0.2                 # Coefficent of friction
    E::FT = 6e6                 # Young's Modulus
end


"""
Constants(::Type{FT}, args...)

A float type FT can be provided as the first argument of any Constants
constructor. A Constants of type FT will be created by passing all other
arguments to the correct constructor. 
"""
Constants(::Type{FT}, args...) where {FT <: AbstractFloat} =
    Constants{FT}(args...)

"""
    Constants(args...)

If a type isn't specified, Constants will be of type Float64 and the correct
constructor will be called with all other arguments.
"""
Constants(args...) = Constants{Float64}(args...)

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
    RT<:Random.AbstractRNG,
    IW<:StructVector{<:InitialStateOutputWriter},
    FW<:StructVector{<:FloeOutputWriter},
    GW<:StructVector{<:GridOutputWriter},
    CW<:StructVector{<:CheckpointOutputWriter},
}
    model::Model{FT, GT, DT}            # Model to simulate
    consts::Constants{FT} = Constants() # Constants used in Simulation
    rng::RT = Xoshiro()                     # Random number generator 
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
    writers::OutputWriters{IW, FW, GW, CW} = OutputWriters()
end

"""
timestep_sim!(sim, tstep, writers)

Run one step of the simulation and write output. 
Inputs:
    sim     <Simulation> simulation to advance
    tstep   <Int> current timestep
    writers <Vector{AbstractOutputWriter}> list of output OutputWriters
Outputs:
    None. Simulation advances by one timestep. 
"""
function timestep_sim!(sim, tstep)
    if !isempty(sim.model.floes)
        max_floe_id = maximum(sim.model.floes.id)
        # Need to lock some operations when multi-threading
        spinlock = Threads.SpinLock()
        # Add ghost floes through periodic boundaries
        n_init_floes = length(sim.model.floes) # number of floes before ghosts
        add_ghosts!(sim.model.floes, sim.model.domain)

        # Output at given timestep
        write_data!(sim, tstep)  # Horribly type unstable
        
        # Collisions
        if sim.collision_settings.collisions_on
            timestep_collisions!(
                sim.model.floes,
                n_init_floes,
                sim.model.domain,
                sim.consts,
                sim.Δt,
                sim.collision_settings,
                spinlock,
            )
        end

        # Remove the ghost floes - only used for collisions
        for i in reverse(n_init_floes+1:length(sim.model.floes))
            StructArrays.foreachfield(
                    field -> deleteat!(field, i),
                    sim.model.floes,
            )
        end
        empty!.(sim.model.floes.ghosts) 

        # Physical processes without ghost floes
        # Effects of ocean and atmosphere on ice and visa versa
        if sim.coupling_settings.coupling_on && mod(tstep, sim.coupling_settings.Δt) == 0
            timestep_coupling!(
                sim.model,
                sim.Δt,
                sim.consts,
                sim.coupling_settings,
            )
        end
        
        # Move and update floes based on collisions and ocean/atmosphere forcing
        timestep_floe_properties!(
            sim.model.floes,
            sim.Δt,
            sim.simp_settings.max_floe_height,
        )

        # TODO: Remove parent ids ?
        # Fracture floes
        if sim.fracture_settings.fractures_on && mod(tstep, sim.fracture_settings.Δt) == 0
            max_floe_id =
                fracture_floes!(
                    sim.model.floes,
                    max_floe_id,
                    sim.rng,
                    sim.fracture_settings,
                    sim.coupling_settings,
                    sim.simp_settings,
                    sim.consts,
                    sim.Δt,
                )
        end
        max_floe_id = 
            simplify_floes!(
                sim.model,
                max_floe_id,
                sim.simp_settings,
                sim.collision_settings,
                sim.coupling_settings,
                sim.Δt,
                sim.consts,
                sim.rng,
            )
    end

    # h0 = real(sqrt.(Complex.((-2Δt * newfloe_Δt) .* hflx)))
    # mean(h0)
    return 
end

"""
    run!(sim, writers)

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
function run!(sim)
    # Start simulation
    sim.verbose && println(string(sim.name, " is running!"))
    tstep = 0
    while tstep <= sim.nΔt
        if sim.verbose && mod(tstep, 50) == 0
            println(tstep, " timesteps")
        end
        # Timestep the simulation forward
        timestep_sim!(sim, tstep)
        tstep+=1
    end
    sim.verbose && println(string(sim.name, " done running!"))
    return
end