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
    model::Model{FT, GT, DT}            # Model to simulate
    consts::Constants{FT}           # Constants used in simulation
    name::String = "sim"            # Simulation name for saving output
    Δd::Int = 1                     # Number of buffer grid cells on each side of floe for monte carlo interpolation
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

function ghosts_on_bounds(_, ghosts, ::NonPeriodicBoundary, _)
    return ghosts
end

function ghosts_on_bounds(primary, ghosts, boundary::PeriodicBoundary, trans_vec)
    if LG.intersects(LG.Polygon(primary.coords), LG.Polygon(boundary.coords))
        # ghosts of existing ghosts and original element
        for elem in vcat(ghosts, primary)
            g = deepcopy(elem)
            g.coords = translate(g.coords, trans_vec)
            g.centroid = g.centroid .+ trans_vec
            push!(ghosts, g)
        end
    end
    return ghosts
end

function find_ghosts(f, domain, ::Type{T} = Float64) where T
    ghosts = Vector{Floe}()
    # Add east/west ghost floes
    Lx = domain.east.val - domain.west.val
    if f.centroid[1] - f.rmax < domain.west.val
        ghosts = ghosts_on_bounds(f, ghosts, domain.west, [Lx, T(0)])
    elseif (f.centroid[1] + f.rmax > domain.east.val)
        ghosts = ghosts_on_bounds(f, ghosts, domain.east, [-Lx, T(0)])
    end
    # Add north/south ghost floes
    Ly =  domain.north.val - domain.south.val
    if (f.centroid[2] - f.rmax < domain.south.val)
        ghosts = ghosts_on_bounds(f, ghosts, domain.south, [T(0), Ly])
    elseif (f.centroid[2] + f.rmax > domain.north.val)
        ghosts = ghosts_on_bounds(f, ghosts, domain.north, [T(0), -Ly])
    end
    println(typeof(ghosts))
    return f, ghosts
end

function add_ghost_floes!(floes, domain, ::Type{T} = Float64) where T
    nfloes = length(floes)
    for i in eachindex(floes)  # uses initial length of floes so we can append to list
        f = floes[i]
        if f.alive == 1
            ghosts = find_ghosts(f, domain, T)
            # If we have ghost floes, update lists
            if !isempty(ghosts)
                inbound_idx = findfirst(g -> (domain.west.val < g.centroid[1] < domain.east.val && domain.south.val < g.centroid[2] < domain.north.val), ghosts)
                if !isnothing(inbound_idx)
                    println(typeof(f))
                    println(typeof(ghosts))
                    temp = ghosts[inbound_idx]
                    ghosts[inbound_idx] = f
                    f = temp
                    #f, ghosts[inbound_idx] = ghosts[inbound_idx], f  # "primary" floe now has centroid within domain
                end
                ghosts.id *= -1  # ghosts have negative index of parent floe
                append!(floes, ghosts)
                append!(f.ghosts, nfloes+1:nfloes+length(ghosts))  # index of ghosts floes saved
                nfloes += length(ghosts)
                floes[i] = f
            end
        end
    end
end

function add_ghost_topography!(domain, ::Type{T} = Float64) where T
    topography = domain.topography
    for i in eachindex(topography)  # uses initial length of topography so we can append to list
        ghosts = find_ghosts(topography[i], domain, T)
        if !isempty(ghosts)
            append!(domain.topography, ghosts)
        end
    end
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
    floe.collision_force[1] += sum(floe.interactions[:, "xforce"])
    floe.collision_force[2] += sum(floe.interactions[:, "yforce"])
    floe.collision_trq += sum(floe.interactions[:, "torque"])

    cforce = floe.collision_force
    ctrq = floe.collision_trq

    if floe.height > 10
        floe.height = 10
    end
    # TODO: Make variable to user input
    if floe.mass < 100
        floe.mass = 1e3
        floe.alive = 0
    end

    while maximum(abs.(cforce)) > floe.mass/(5Δt)
        cforce = cforce ./ 10
        ctrq = ctrq ./ 10
        # TODO: check floe interactions
    end
    h = floe.height
    # Update floe based on thermodynamic growth
    Δh = floe.hflx * Δt/h
    Δh = 0 # for collision testing
    hfrac = (h-Δh)/h
    floe.mass *= hfrac
    floe.moment *= hfrac
    floe.height -= Δh
    h = floe.height

    # Update ice coordinates with velocities and rotation
    Δx = 1.5Δt*floe.u - 0.5Δt*floe.p_dxdt
    Δy = 1.5Δt*floe.v - 0.5Δt*floe.p_dydt
    floe.centroid .+= [Δx, Δy]
    floe.coords = translate(floe.coords, [Δx, Δy])
    floe.p_dxdt = floe.u
    floe.p_dydt = floe.v

    Δα = 1.5Δt*floe.ξ - 0.5Δt*floe.p_dαdt
    floe.α += Δα
    floe.p_dαdt = floe.ξ
    floe.coords = [map(p -> [cos(Δα)*p[1] - sin(Δα)*p[2],
                              sin(Δα)*p[1] + cos(Δα)p[2]], floe.coords[1])]

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

"""
    timestep_atm!(m)

Update model's ocean and heat flux from atmosphere effects. 
Input:
        m <Model>
Outputs: 
        None. The ocean stress fields are updated from atmos stress.
        Heatflux is also updated with ocean and atmos temperatures. 
"""
function timestep_atm!(m, c)
    Δu_AO = m.atmos.u .- m.ocean.u
    Δv_AO = m.atmos.v .- m.ocean.v
    m.ocean.taux .= c.ρa  *c.Cd_ao * sqrt.(Δu_AO.^2 + Δv_OI.^2) .* Δu_AO
    m.ocean.tauy .= c.ρa * c.Cd_ao * sqrt.(Δu_AO.^2 + Δv_OI.^2) .* Δv_AO
    m.hflx .= c.k/(c.ρi*c.L) .* (atmos.temp .- ocean.temp)
end

function timestep_sim!(sim, tstep, ::Type{T} = Float64) where T
    m = sim.model
    n_init_floes = length(m.floes) # number of floes before ghost floes
    m.ocean.si_area .= zeros(T, 1)
    remove = zeros(Int, n_init_floes)
    transfer = zeros(Int, n_init_floes)
    add_ghost_floes!(m.floes, m.domain)
    for i in eachindex(m.floes) # floe-floe collisions for floes i and j where i<j
        ifloe = m.floes[i]
        ifloe.collision_force = zeros(T, 1, 2)
        ifloe.collision_trq = T(0.0)
        ifloe.interactions = ifloe.interactions[1:1, :]
        if sim.COLLISION
            for j in i+1:length(floes)
                if sum((ifloe.centroid .- m.floes[j].centroid).^2) < (ifloe.rmax + m.floes[j].rmax)^2
                    ikill, itransfer = floe_floe_interaction!(ifloe, i, m.floes[j], j, n_init_floes, sim.consts, sim.Δt)
                    if ikill != 0 || itransfer != 0
                        remove[i] = ikill
                        transfer[i] = itransfer
                    end
                end
            end
        end
        floe_domain_interaction!(ifloe, m.domain, sim.consts, sim.Δt)
        m.floes[i] = ifloe
    end
    for i in eachindex(m.floes)  # Update floes not directly calculated above where i>j
        if i <= n_init_floes && remove[i] > 0 && remove[i] != i
            transfer[remove[i]] = i
        end
        ij_inters = m.floes[i].interactions
        if size(ij_inters, 1) > 1 # First row is place holder for initialization and clearing
            for inter_idx in axes(ij_inters, 1)  # Loop over each interaction with Floe i
                j = ij_inters[inter_idx, "floeidx"]  # Index of floe to update in model floe list
                if j <= length(m.floes) && j > i
                    jidx = Int(j)
                    jfloe = m.floes[jidx]
                    jfloe.interactions = [jfloe.interactions; i -ij_inters[inter_idx, "xforce"] -ij_inters[inter_idx, "yforce"] #=
                                        =# ij_inters[inter_idx, "xpoint"] ij_inters[inter_idx, "ypoint"] T(0.0) ij_inters[inter_idx, "overlap"]]
                    jfloe.overarea += ij_inters[inter_idx, "overlap"]
                    m.floes[jidx] = jfloe
                end
            end
        end
    end
    for i in reverse(eachindex(m.floes)) # TODO: Move torque calculations to earlier and then we only need to do this for the first points
        ifloe = m.floes[i]
        calc_torque!(ifloe)
        for g in ifloe.ghosts
            gfloe = m.floes[g]
            ifloe.collision_force[1] += sum(gfloe.interactions[:, "xforce"])
            ifloe.collision_force[2] += sum(gfloe.interactions[:, "yforce"])
            ifloe.collision_trq += sum(gfloe.interactions[:, "torque"])

        end
        m.floes[i] = ifloe
    end
    m.floes = m.floes[1:n_init_floes] # remove the ghost floes
    for i in 1:n_init_floes
        ifloe = m.floes[i]
        if mod(tstep, sim.Δtocn) == 0
            floe_OA_forcings!(ifloe, m, sim.consts, sim.Δd)
        end
        timestep_floe!(ifloe, sim.Δt)
        m.floes[i] = ifloe
    end
    
    remove_idx = findall(f -> f.alive == 0, m.floes)
    for idx in remove_idx
        StructArrays.foreachfield(f -> deleteat!(f, idx), m.floes)
    end
    m.ocean.hflx .= sim.consts.k/(sim.consts.ρi*sim.consts.L) .* (m.atmos.temp .- m.ocean.temp)
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
    # Output setup
    Δtout_lst = Int[]
    if !isempty(writers)
        write_domain!(sim.model.domain, sim.name)
    end
    for w in writers
        setup_output_file!(w, sim.nΔt, sim.name, T)
        push!(Δtout_lst, w.Δtout)
    end
    println(string(sim.name ," running!"))

    # Initialize floe IDs
    for i in eachindex(sim.model.floes)
        f = sim.model.floes[i]
        f.id = i
        sim.model.floes[i] = f
    end
    
    # Add topography elements crossing through periodic boundaries
    add_ghost_topography!(sim.model.domain)

    tstep = 0
    while tstep <= sim.nΔt
        # Write data
        widx = findall(Δtout-> mod(tstep, Δtout) == 0, Δtout_lst)
        if length(widx) > 0
            println(tstep, " timesteps completed")
            for idx in widx
                write_data!(writers[idx], tstep, sim.model, sim.name)
            end
        end
        # Timestep the simulation forward
        timestep_sim!(sim, tstep, T)
        tstep+=1
    end
    println(string(sim.name ," done running!"))
end