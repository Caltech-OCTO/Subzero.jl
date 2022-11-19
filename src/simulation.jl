"""
Structs and functions to create and run a Subzero simulation
"""

"""
    Simulation{FT<:AbstractFloat, DT<:AbstractDomain{FT}}

Simulation which holds a model and parameters needed for running the simulation. Simulation requires a model, a coarse grid, a coarse grid data struct, and a figure. The figure can be initialized using setup_plot. The rest of the simulation values are optional. These fields and their default values are as follows:  the size of a timestep in seconds Δt (10), the total number of timesteps in the simulation nΔt (7500), the output frequency of floe and data on the coarse grid in timesteps nΔtout (150),  timesteps between saving images Δtpics (150),  timesteps between floe simplicaiton  Δtsimp (20), timesteps betwen thermodynamic floe creation Δtpack (500), timesteps between updating ocean forcing  Δtocn (10). There are also flags that control simulation behavior. These flags are AVERAGE (average coarse grid data in time), COLLISION (enable floe collisions), CORNERS (floe corners can break), FRACTURES (floes can fracture), KEEPMIN (small floes don't dissolve), PACKING (floe packing enabled), RAFTING (floe rafting enabled), RIDGING (floe ridging enabled), and WELDING (floe welding enabled). All are false by default.
"""
@kwdef struct Simulation{FT<:AbstractFloat, DT<:AbstractDomain{FT}}
    # Objects ------------------------------------------------------------------
    model::Model{FT, DT}            # Model to simulate
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
        floe with updated fields
"""
function timestep_floe!(floe, Δt)
    floe.collision_force[1] += sum(floe.interactions[:, "xforce"])
    floe.collision_force[2] += sum(floe.interactions[:, "yforce"])
    floe.collision_torque += sum(floe.interactions[:, "torque"])

    cforce = floe.collision_force
    ctorque = floe.collision_torque

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
        ctorque = ctorque ./ 10
        # TODO: check floe interactions
    end
    h = floe.height
    # Update floe based on thermodynamic growth
    Δh = floe.hflx * Δt/h
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

    dξdt = (floe.torqueOA + ctorque)/floe.moment
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

function timestep_atm!(m)
    c = m.constants
    Δu_AO = m.wind.u .- m.ocean.u
    Δv_AO = m.wind.v .- m.ocean.v
    m.ocean.taux .= c.ρa  *c.Cd_ao * sqrt.(Δu_AO.^2 + Δv_OI.^2) .* Δu_AO
    m.ocean.tauy .= c.ρa * c.Cd_ao * sqrt.(Δu_AO.^2 + Δv_OI.^2) .* Δv_AO
    m.hflx .= c.k/(c.ρi*c.L) .* (wind.temp .- ocean.temp)
end


function floe_boundary_interaction!(floe, boundary_coords, ::OpenBC, _)
    floe_poly = LG.Polygon(floe.coords)
    bounds_poly = LG.Polygon(boundary_coords)
    # Check if the floe and boundary actually overlap
    if LG.intersects(floe_poly, bounds_poly)
        floe.alive = 0
    end
    return
end

function floe_boundary_interaction!(floe, boundary_coords, ::CollisionBC,
consts, t::Type{T} = Float64) where T
    # Check if the floe and boundary actually overlap
    floe_poly = LG.Polygon(floe.coords)
    bounds_poly = LG.Polygon(boundary_coords)
    if LG.intersects(floe_poly, bounds_poly)
        inter_floe = LG.intersection(floe_poly, bounds_poly)
        inter_regions = LG.getGeometries(inter_floe)
        region_areas = [LG.area(poly) for poly in inter_regions]::Vector{Float64}
        if maximum(region_areas)/floe.area > 0.75  # Regions overlap too much
            return zeros(T, 1, 2), zeros(T, 1, 2), fill(T(Inf), 1)
        end
        # Constant needed for force calculations
        force_factor = consts.E * floe.height / sqrt(floe.area)
        # Total collision forces
        forces, pcontacts, overlaps =  calc_collision_forces(floe.coords, boundary_coords,
                                        inter_regions, region_areas, force_factor, consts, T)
        # Add this back in once we figure out Ly, Lx --> should we just use poly2 to check if ON
        if sum(abs.(forces)) != 0
        #   forces[pcontacts[:, 2] .== Ly, 1] .= T(0.0)
        #   forces[pcontacts[:, 1] .== Lx, 2] .= T(0.0)
            nint = size(forces, 1)
            floe.interactions = [floe.interactions; fill(Inf, nint) forces pcontacts zeros(nint) overlaps']
            floe.overarea += sum(overlaps)
        end
    end
    return
end

function floe_domain_interaction!(floe, domain::DT, consts, t::Type{T} = Float64) where {DT<:RectangleDomain, T}
    centroid = floe.centroid
    rmax = floe.rmax
    eastval = domain.east.val
    westval = domain.west.val
    northval = domain.north.val
    southval = domain.south.val

    # unless periodic
    if centroid[1] > eastval || centroid[1] < westval || centroid[2] > northval || centroid[2] < southval
        floe.alive = 0
        return
    else

    end
   
    if centroid[1] + rmax > eastval
        boundary_coords = [[[eastval, northval], [eastval, southval],
                            [eastval + floe.rmax, southval],
                            [eastval + rmax, northval], [eastval, northval]]]
        floe_boundary_interaction!(floe, boundary_coords, domain.east.bc, consts)
    end
    if centroid[1] - rmax < westval
        boundary_coords = [[[westval, northval], [westval, southval],
                            [westval - floe.rmax, southval],
                            [westval - rmax, northval], [westval, northval]]]
        floe_boundary_interaction!(floe, boundary_coords, domain.west.bc, consts)
    end
    if centroid[2] + rmax > northval
        boundary_coords = [[[westval, northval], [westval, northval + rmax],
                            [eastval, northval + floe.rmax],
                            [eastval, northval], [westval, northval]]]
        floe_boundary_interaction!(floe, boundary_coords, domain.north.bc, consts)
    end
    if centroid[2] - rmax < southval
        boundary_coords = [[[westval, southval], [westval, southval - rmax],
                            [eastval, southval - floe.rmax],
                            [eastval, southval], [westval, southval]]]
        floe_boundary_interaction!(floe, boundary_coords, domain.north.bc, consts)
    end
    return
end

function floe_floe_interaction!(ifloe, i, jfloe, j, nfloes, consts, t::Type{T} = Float64) where T
    remove = T(0.0)
    transfer = T(0.0)
    ifloe_poly = LG.Polygon(ifloe.coords)
    jfloe_poly = LG.Polygon(jfloe.coords)
    if LG.intersects(ifloe_poly, jfloe_poly)  # Interactions
        inter_floe = LG.intersection(ifloe_poly, jfloe_poly)
        inter_regions = LG.getGeometries(inter_floe)
        region_areas = [LG.area(poly) for poly in inter_regions]::Vector{Float64}
        total_area = sum(region_areas)
        # Floes overlap too much - remove floe or transfer floe mass to other floe
        if total_area/ifloe.area > 0.55
            if i <= nfloes
                remove = i
                transfer = j
            elseif j <= nfloes
                remove = j  # Will transfer mass to ifloe 
            end
        elseif total_area/jfloe.area > 0.55
            if j <= nfloes
                remove = j  # Will transfer mass to ifloe
            end
        else
            # Calculate force force
            ih = ifloe.height
            ir = sqrt(ifloe.area)
            jh = jfloe.height
            jr = sqrt(jfloe.area)
            force_factor = if ir>1e5 || jr>1e5
                consts.E*min(ih, jh)/min(ir, jr)
            else
                consts.E*(ih*jh)/(ih*jr+jh*ir)
            end
            # Calculate forces, force points, and overlap areas
            forces, pcontacts, overlaps = calc_collision_forces(ifloe.coords, jfloe.coords,
                                        inter_regions, region_areas, force_factor, consts, T)
            if sum(abs.(forces)) != 0
                nint = size(forces, 1)
                ifloe.interactions = [ifloe.interactions; fill(j, nint) forces pcontacts zeros(nint) overlaps']
                ifloe.overarea += sum(overlaps)
            end
        end
    end
    return remove, transfer
end

"""
    domain_coords(domain::RectangleDomain)
Inputs:
        domain<RectangleDomain>
Output:
        RingVec coordinates for edges of rectangular domain based off of boundary values
"""
function cell_coords(xmin, xmax, ymin, ymax)
    return [[[xmin, ymax], [xmin, ymin],
             [xmax, ymin], [xmax, ymax],
             [xmin, ymax]]]
end

function run!(sim, writers, t::Type{T} = Float64) where T
    Δtout_lst = Int[]
    for w in writers
        setup_output_file!(w, sim.nΔt, T)
        push!(Δtout_lst, w.Δtout)
    end

    println("Model running!")
    m = sim.model
    plt = setup_plot(m)
    tstep = 0
    plot_sim(m, plt, tstep)
    while tstep <= sim.nΔt
        widx = findall(Δtout-> mod(tstep, Δtout) == 0, Δtout_lst)
        if length(widx) > 0
            println(tstep, " timesteps completed")
            for idx in widx
                write_data!(writers[idx], tstep, m)
            end
            plot_sim(m, plt, tstep)
        end
        m.ocean.si_frac .= zeros(T, 1)
        nfloes = length(m.floes) # number of floes before ghost floes
        remove = zeros(Int, nfloes)
        transfer = zeros(Int, nfloes)
        for i in 1:nfloes # floe-floe collisions for floes i and j where i<j
            ifloe = m.floes[i]
            ifloe.collision_force = zeros(T, 1, 2)
            ifloe.collision_torque = T(0.0)
            ifloe.interactions = ifloe.interactions[1:1, :]
            if sim.COLLISION
                for j in i+1:nfloes
                    if sum((ifloe.centroid .- m.floes[j].centroid).^2) < (ifloe.rmax + m.floes[j].rmax)^2
                        ikill, itransfer = floe_floe_interaction!(ifloe, i, m.floes[j], j, nfloes, m.consts)
                        remove[i] = ikill
                        transfer[i] = itransfer
                    end
                end
            end
            m.floes[i] = ifloe
        end
        for i in 1:length(m.floes)  # Update floes not directly calculated above where i>j
            if i <= nfloes && remove[i] > 0 && remove[i] != i
                transfer[remove[i]] = j
            end
            ij_inters = m.floes[i].interactions
            if size(ij_inters, 1) > 1
                for j in ij_inters[:, "floeidx"]
                    if j <= length(m.floes) && j > i
                        jidx = Int(j)
                        jfloe = m.floes[jidx]
                        jfloe.interactions = [jfloe.interactions; i -ij_inters[jidx, "xforce"] -ij_inters[jidx, "yforce"] #=
                                           =# ij_inters[jidx, "xpoint"] ij_inters[jidx, "ypoint"] T(0.0) ij_inters[jidx, "overlap"]]
                        jfloe.overarea += ij_inters[jidx, "overlap"]
                        m.floes[jidx] = jfloe
                    end
                end
            end
        end
        for i in 1:nfloes
            ifloe = m.floes[i]
            floe_domain_interaction!(ifloe, m.domain, m.consts)
            floe_OA_forcings!(ifloe, m)
            #calc_torque!(ifloe)
            timestep_floe!(ifloe, sim.Δt)
            m.floes[i] = ifloe
        end
        

        remove_idx = findall(f -> f.alive == 0, m.floes)
        for idx in remove_idx
            StructArrays.foreachfield(f -> deleteat!(f, idx), m.floes)
        end
        tstep+=1
    end

    # h0 = real(sqrt.(Complex.((-2Δt * newfloe_Δt) .* hflx)))
    # mean(h0)
    plot_sim(m, plt, tstep)
    println("Model done running!")
end