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

Simulation which holds a model and parameters needed for running the simulation. Simulation requires a model, a coarse grid, a coarse grid data struct, and a figure.
    The figure can be initialized using setup_plot. The rest of the simulation values are optional. These fields and their default values are as follows:
    the size of a timestep in seconds Δt (10), the total number of timesteps in the simulation nΔt (7500), the output frequency of floe and data on the coarse grid in
    timesteps nΔtout (150),  timesteps between saving images Δtpics (150),  timesteps between floe simplicaiton  Δtsimp (20), timesteps betwen thermodynamic floe creation Δtpack (500),
    timesteps between updating ocean forcing  Δtocn (10). There are also flags that control simulation behavior. These flags are AVERAGE (average coarse grid data in time),
    COLLISION (enable floe collisions), CORNERS (floe corners can break), FRACTURES (floes can fracture), KEEPMIN (small floes don't dissolve), PACKING (floe packing enabled),
    RAFTING (floe rafting enabled), RIDGING (floe ridging enabled), and WELDING (floe welding enabled). All are false by default.
"""
@kwdef struct Simulation{FT<:AbstractFloat, GT<:AbstractGrid, DT<:Domain}
    model::Model{FT, GT, DT}        # Model to simulate
    consts::Constants{FT}           # Constants used in Simulation
    Δd::Int = 1                     # Number of buffer grid cells on each side of floe for monte carlo interpolation
    verbose::Bool = false           # String output printed during run
    # Timesteps ----------------------------------------------------------------
    Δt::Int = 10                    # Simulation timestep (seconds)
    nΔt::Int = 7500                 # Total timesteps simulation runs for
    Δtsimp::Int = 20                # Timesteps between floe simplification
    Δtpack::Int = 500               # Timesteps between thermodynamic floe 
                                    # creation
    Δtocn::Int = 10                 # Timesteps between updating ocean forces
    # Flags --------------------------------------------------------------------
    COLLISION::Bool = false         # If true, collisions are enabled for floes
    CORNERS::Bool = false           # If true, corners of floes can break
    FRACTURES::Bool = false         # If true, fracturing of floes is enabled
    KEEPMIN::Bool = false           # If true, retain small floes that would 
                                    # normally "dissolve"
    PACKING::Bool = false           # If true, floe packing is enabled
    RAFTING::Bool = false           # If true, floe rafting is enabled
    RIDGING::Bool = false           # If true, floe ridging is enabled
    WELDING::Bool = false           # If true, floe welding is enabled
end

"""
    timestep_floe(floe)

Update floe position and velocities using second-order time stepping with tendencies calculated at previous timesteps.
Input:
        floe <Floe>
Output:
        None. Floe's fields are updated with new position and speed
"""
function timestep_floe!(floe, Δt)
    cforce = floe.collision_force
    ctrq = floe.collision_trq

    if floe.height > 10
        floe.height = 10
    end
    # TODO: Make variable to user input
    if floe.mass < 100
        floe.mass = 1e3
        floe.alive = false
    end

    while maximum(abs.(cforce)) > floe.mass/(5Δt)
        cforce = cforce ./ 10
        ctrq = ctrq ./ 10
    end
    
    h = floe.height
    # Update floe based on thermodynamic growth
    Δh = floe.hflx_factor / h
    #Δh = 0 # for collision testing
    hfrac = (h + Δh) / h
    floe.mass *= hfrac
    floe.moment *= hfrac
    floe.height -= Δh
    h = floe.height

    # Update ice coordinates with velocities and rotation
    Δx = 1.5Δt*floe.u - 0.5Δt*floe.p_dxdt
    Δy = 1.5Δt*floe.v - 0.5Δt*floe.p_dydt
    Δα = 1.5Δt*floe.ξ - 0.5Δt*floe.p_dαdt
    floe.α += Δα

    coords0 = translate(floe.coords, -floe.centroid)
    coords0 = [map(p -> [cos(Δα)*p[1] - sin(Δα)*p[2],
                         sin(Δα)*p[1] + cos(Δα)p[2]], coords0[1])]
    floe.centroid .+= [Δx, Δy]
    floe.coords = translate(coords0, floe.centroid)

    floe.p_dxdt = floe.u
    floe.p_dydt = floe.v
    floe.p_dαdt = floe.ξ

    # Update ice velocities with forces and torques
    dudt = (floe.fxOA + cforce[1])/floe.mass
    dvdt = (floe.fyOA + cforce[2])/floe.mass
    frac = if abs(Δt*dudt) > (h/2) && abs(Δt*dvdt) > (h/2)
        frac1 = (sign(dudt)*h/2Δt)/dudt
        frac2 = (sign(dvdt)*h/2Δt)/dvdt
        min(frac1, frac2)
    elseif abs(Δt*dudt) > (h/2) && abs(Δt*dvdt) < (h/2)
        (sign(dudt)*h/2Δt)/dudt
    elseif abs(Δt*dudt) < (h/2) && abs(Δt*dvdt) > (h/2)
        (sign(dvdt)*h/2Δt)/dvdt
    else
        1
    end

    dudt = frac*dudt
    dvdt = frac*dvdt
    floe.u += 1.5Δt*dudt-0.5Δt*floe.p_dudt
    floe.v += 1.5Δt*dvdt-0.5Δt*floe.p_dvdt
    floe.p_dudt = dudt
    floe.p_dvdt = dvdt

    dξdt = (floe.trqOA + ctrq)/floe.moment
    dξdt = frac*dξdt
    ξ = floe.ξ + 1.5Δt*dξdt-0.5Δt*floe.p_dξdt
    if abs(ξ) > 1e-5
        ξ = sign(ξ) * 1e-5
    end
    floe.ξ = ξ
    floe.p_dξdt = dξdt
    return
    # TODO: Floe strain - Calc_trajectory lines 216-288
    # TODO: Floe stress - Calc_trajectory lines 9-21
end

function timestep_sim!(sim, tstep, writers, ::Type{T} = Float64) where T
    m = sim.model
    # Add ghosts
    n_init_floes = length(m.floes) # number of floes before ghost floes
    add_ghosts!(m.floes, m.domain)

    # Output at given timestep
    for w in writers
        if tstep == 0 || (hasfield(typeof(w), :Δtout) && mod(tstep, w.Δtout) == 0)
            write_data!(w, tstep, sim)
        end
    end

    # Collisions
    remove = zeros(Int, n_init_floes)
    transfer = zeros(Int, n_init_floes)
    if sim.COLLISION
        remove, transfer = timestep_collisions!(m.floes, n_init_floes, m.domain, remove, transfer, sim.consts, sim.Δt, T)
    end

    m.floes = m.floes[1:n_init_floes] # remove the ghost floes
    empty!.(m.floes.ghosts) 
    for i in 1:n_init_floes
        ifloe = m.floes[i]
        #if mod(tstep-1, sim.Δtocn) == 0
        floe_OA_forcings!(ifloe, m, sim.consts, sim.Δd)
        #end
        timestep_floe!(ifloe, sim.Δt)
        m.floes[i] = ifloe
    end

    # Timestep ocean
    timestep_ocean!(m, sim.consts, sim.Δt)
    
    remove_idx = findall(f -> !f.alive, m.floes)
    for idx in remove_idx
        StructArrays.foreachfield(f -> deleteat!(f, idx), m.floes)
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
        None. The simulation will be run and outputs will be saved in the output folder. 
"""
function run!(sim, writers, ::Type{T} = Float64) where T
    # Start simulation
    sim.verbose && println("Simulation running!")
    tstep = 0
    while tstep <= sim.nΔt
        if sim.verbose && mod(tstep, 50) == 0
            println(tstep, " timesteps")
        end
        # Timestep the simulation forward
        timestep_sim!(sim, tstep, writers, T)
        tstep+=1
    end
    sim.verbose && println("Simulation done running!")
end