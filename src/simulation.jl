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
    floe_grid_bounds(g, p, rmax)

Finds the bounding indices of a grid line such that a point plus and minus a maximum radius are within those grid lines.
Inputs:
        g <Vector{Real}> grid lines
        p <Real> point value
        rmax <Real> radial buffer around p
Outputs:
        min_idx: index of g such that the value of (p - rmax) is between 
                 indices g[min_idx] and g[min_idx + 1]. If p is less than all of the values in g, this will be 0.
        max_idx: index of g such that the value of (p + rmax) is between 
                 indices g[max_idx - 1] and g[man_idx]. If p is greater than all of the values in g, this will be the length of g.
Note: If radius is negative this will switch the minimum and maximum indices.
"""
function floe_grid_bounds(g, p, rmax)
    Δ = g[2] - g[1]
    pmin = p - rmax > g[1] ? p - rmax : g[1]
    pmax = p + rmax < g[end] ? p + rmax : g[end]
    min_idx = findmin(abs.(g .- pmin))[2]
    min_val = g[min_idx]
    if pmin < min_val
        min_idx -= 1
        min_val -= Δ
    end
    max_idx = Int(cld(pmax-pmin, Δ)) + min_idx
    if pmax > g[max_idx]
        max_idx += 1
    end
    return min_idx, max_idx
end

"""
    cell_area_ratio(cell_coords, floe_coords)

Calculates the percentage of a grid square filled with a given floe.
Inputs:
        cell_poly <LibGEOS.Polygon>
        floe_poly <LibGEOS.Polygon>
Outputs:
        cell area ratio <Float> ratio of cell area filled with given floe
"""
function cell_area_ratio(cell_poly, floe_poly)
    floe_in_cell = LG.intersection(floe_poly, cell_poly)
    return LG.area(floe_in_cell)/LG.area(cell_poly)
end

"""
    floe_area_ratio(floe, xg, yg)

Calculates the cell area ratio of grid squares surrounding given floe and the indicies of those grid squares within the grid defined by gridlines xg and yg.
Inputs:
        floe    <Floe>
        xg      <Vector{Float}> x grid lines
        yg      <Vector{Float}> y grid lines
Outputs:
        area_ratios <Vector{Float}> vector of area ratios for each grid cell 
                                    within floe grid bounds
        xidx <Vector{Int}> x indices of grid cells - order matches area_ratios
        yidx <Vector{Int}> y indices of grid cells - order matches area_ratios
        idx <Vector{(Int, Int)}> cartesian point defining one grid cell in grid
                                 - order matches area_ratios
"""
function floe_area_ratio(floe, xg, yg, t::Type{T} = Float64) where T
    floe_poly = LG.Polygon(translate(floe.αcoords, floe.centroid))
    xmin_idx, xmax_idx = floe_grid_bounds(xg, floe.centroid[1], floe.rmax)
    ymin_idx, ymax_idx = floe_grid_bounds(yg, floe.centroid[2], floe.rmax)
    nx = xmax_idx - xmin_idx
    ny = ymax_idx - ymin_idx
    area_ratios = T[]
    xidx = Int[]
    yidx = Int[]
    for i = xmin_idx:(xmax_idx-1)
        for j = ymin_idx:(ymax_idx-1)
            cell_poly = LG.Polygon(cell_coords(xg[i], xg[i+1], yg[j], yg[j+1]))
            ratio = cell_area_ratio(cell_poly, floe_poly)
            if ratio > 0.0
                push!(area_ratios, ratio)
                push!(xidx, i)
                push!(yidx, j)
            end
        end
    end
    idx = CartesianIndex.(Tuple.(eachrow(hcat(xidx,yidx))))
    return area_ratios, xidx, yidx, idx
end

"""
    calc_OA_forcings!(m, i)

Calculate the effects on the ocean and atmpshere on floe i within the given model
and the effects of the ice floe on the ocean.

Inputs:
        m <Model> given model
        i <Int> floe i within the model's floes field
Outputs:
        None. Both floe and ocean fields are updated in-place. 
"""
function calc_OA_forcings!(m, i)
    floe = m.floes[i]
    Δx = m.grid.xg[2] - m.grid.xg[1]
    Δy = m.grid.yg[2] - m.grid.yg[1]

    # Grid squares under ice floe and ice area per cell
    ma_ratio = floe.mass/floe.area
    area_ratios, xidx, yidx, idx = floe_area_ratio(floe, m.grid.xg, m.grid.yg)
    areas = area_ratios * (Δx * Δy)

    # Floe heatflux
    floe.hflx = mean(m.hflx[idx])

    # Ice velocity within each grid square
    lx = m.grid.xc[xidx] .- floe.centroid[1]
    ly = m.grid.yc[yidx] .- floe.centroid[2]
    uice = floe.u .- ly*floe.ξ
    vice = floe.v .- lx*floe.ξ

    # Force on ice from atmopshere
    uatm = m.wind.u[idx]
    vatm = m.wind.u[idx]
    fx_atm = (m.ρa * m.CIA * sqrt.(uatm.^2 + vatm.^2) .* uatm) .* areas
    fy_atm = (m.ρa * m.CIA * sqrt.(uatm.^2 + vatm.^2) .* vatm) .* areas

    # Force on ice from pressure gradient
    fx_pressure∇ = -ma_ratio * m.coriolis .* m.ocean.v[idx] .* areas
    fy_pressure∇ = ma_ratio * m.coriolis .* m.ocean.u[idx] .* areas

    # Force on ice from ocean
    Δu_OI = m.ocean.u[idx] .- uice
    Δv_OI = m.ocean.v[idx] .- vice
    τx_ocn = m.ρo*m.CIO*sqrt.(Δu_OI.^2 + Δv_OI.^2) .* (cos(m.turnθ) .* Δu_OI .+ sin(m.turnθ) * Δv_OI)
    τy_ocn = m.ρo*m.CIO*sqrt.(Δu_OI.^2 + Δv_OI.^2) .* (-sin(m.turnθ) .* Δu_OI .+ cos(m.turnθ) * Δv_OI)
    fx_ocn = τx_ocn .* areas
    fy_ocn = τy_ocn .* areas

    # Sum above forces and find torque
    fx = fx_atm .+ fx_pressure∇ .+ fx_ocn
    fy = fy_atm .+ fy_pressure∇ .+ fy_ocn
    trq = lx.*fy .- ly.*fx

    # Add coriolis force to total foces
    fx .+= ma_ratio * m.coriolis * floe.v * areas
    fy .-= ma_ratio * m.coriolis * floe.u * areas

    # Sum forces on ice floe
    floe.fxOA = sum(fx)
    floe.fyOA = sum(fy)
    floe.torqueOA = sum(trq)
    m.floes[i] = floe

    # Update ocean stress fields with ice on ocean stress
    m.ocean.τx[idx] .= m.ocean.τx[idx].*(1 .- area_ratios) .- τx_ocn.*area_ratios
    m.ocean.τy[idx] .= m.ocean.τy[idx].*(1 .- area_ratios) .- τy_ocn.*area_ratios

    # Update sea-ice fraction
    m.ocean.si_frac[idx] .+= area_ratios
    return
end

"""
    update_floe!(floe)

Update floe position and velocities using second-order time stepping with tendencies calculated at previous timesteps.
Input:
        floe <Floe>
Output:
        floe with updated fields
"""
function timestep_floe(floe, Δt)
    collision_force = [0.0 0.0]
    collision_torque = 0.0

    if floe.height > 10
        floe.height = 10
    end
    # TODO: Make variable to user input
    if floe.mass < 100
        floe.mass = 1e3
        floe.alive = false
    end

    while maximum(abs.(collision_force)) > floe.mass/(5Δt)
        collision_force /= 10
        collision_torque /= 10
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
    floe.p_dxdt = floe.u
    floe.p_dydt = floe.v

    floe.α += 1.5Δt*floe.ξ - 0.5Δt*floe.p_dαdt
    floe.p_dξdt = floe.ξ
    α = floe.α
    floe.αcoords = [map(p -> [cos(α)*p[1] - sin(α)*p[2],
                              sin(α)*p[1] + cos(α)p[2]], floe.coords[1])]

    # Update ice velocities with forces and torques
    dudt = (floe.fxOA + collision_force[1])/floe.mass
    dvdt = (floe.fyOA + collision_force[2])/floe.mass
    
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

    dξdt = (floe.torqueOA + collision_torque)/floe.moment
    dξdt = frac*dξdt
    ξ = floe.ξ + 1.5Δt*dξdt-0.5Δt*floe.p_dξdt
    if abs(ξ) > 1e-5
        ξ = sign(ξ) * 1e-5
    end
    floe.ξ = ξ
    floe.p_dξdt = dξdt

    return floe

    # TODO: Floe strain - Calc_trajectory lines 216-288
    # TODO: Floe stress - Calc_trajectory lines 9-21
end

# function calc_interations(grid, floe_arr, topography_arr, coarse_nx, coarse_ny, PERIODIC)
#     live_floes = filter(floe->floe.alive, floe_arr)
#     topography_poly = LG.MultiPolygon([translate(poly.coords, poly.centroid) for
#                                        poly in topography_arr])
#     # Periodic section
#     # Multiple floe interactions
# end

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

function run!(sim)
    println("Model running!")
    plt = setup_plot(sim.model)
    t = 1
    plot_sim(sim.model, plt, t)
    while t < sim.nΔt
        fill!(sim.model.ocean.si_frac, 0.0)
        for i in eachindex(sim.model.floes)
            calc_OA_forcings!(sim.model, i)
            new_floe = timestep_floe(sim.model.floes[i], sim.Δt)
            sim.model.floes[i] = new_floe
        end
        println(t)
        t+=1
    end
    plot_sim(sim.model, plt, t)
    println("Model done running!")
end