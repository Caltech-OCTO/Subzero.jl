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
function floe_area_ratio(floe, xg, yg)
    floe_poly = LG.Polygon(translate(floe.αcoords, floe.centroid))
    xmin_idx, xmax_idx = floe_grid_bounds(xg, floe.centroid[1], floe.rmax)
    ymin_idx, ymax_idx = floe_grid_bounds(yg, floe.centroid[2], floe.rmax)
    nx = xmax_idx - xmin_idx
    ny = ymax_idx - ymin_idx
    area_ratios = zeros(nx*ny)
    xidx = zeros(Int, nx*ny)
    yidx = zeros(Int, ny*ny)
    cnt = 1
    for i = xmin_idx:(xmax_idx-1)
        for j = ymin_idx:(ymax_idx-1)
            cell_poly = LG.Polygon(cell_coords(xg[i], xg[i+1], yg[j], yg[j+1]))
            area_ratios[cnt] = cell_area_ratio(cell_poly, floe_poly)
            xidx[cnt] = i
            yidx[cnt] = j
            cnt+=1
        end
    end
    idx = CartesianIndex.(Tuple.(eachrow(hcat(xidx,yidx))))
    return area_ratios, xidx, yidx, idx
end

function calc_OA_forcings!(m, i, Δt, Δtocn)
    floe = m.floe[i]
    ocn = m.ocean
    atm = m.wind
    (Δx, Δy) = m.grid.dims
    area_ratios, xidx, yidx, idx = floe_area_ratio(floe, m.grid.xg, m.grid.yg)

    # Ice stress on ocean
    lx = m.grid.xc[xidx] .- floe.centroid[1]
    ly = m.grid.yc[yidx] .- floe.centroid[2]
    ui = floe.u .- ly*floe.ξ
    vi = floe.v .- lx*floe.ξ
    τxIO = m.ρo*m.Cio*(ui .- m.ocean.u[idx]).*abs.(ui .- m.ocean.u[idx])
    τyIO = m.ρo*m.Cio*(vi .- m.ocean.v[idx]).*abs.(vi .- m.ocean.v[idx])
    ocn.τx[idx] .= ocn.τx[idx].*(1 .- area_ratios) .+ τxIO.*area_ratios
    ocn.τy[idx] .= ocn.τy[idx].*(1 .- area_ratios) .+ τyIO.*area_ratios

    # Forces and torques on floes
    Fxx = (-τxIO) + m.ρa * m.CIA * (atm.u[idx_floe] .- ui) .*
          abs.(atm.u .- ui) .* area_ratios*(Δx^2)
    Fyy = (-τyIO) + m.ρa * m.CIA * (atm.v[idx_floe] .- vi) .*
          abs.(atm.v .- vi) .* area_ratios*(Δy^2)
    floe.fxOA = sum(Fxx)
    floe.fyOA = sum(Fyy)
    floe.torqueOA = sum(lx.*Fyy .- ly.*Fxx)
end

"""
    update_floe!(floe)

Update floe position and velocities using second-order time stepping with tendencies calculated at previous timesteps.
Input:
        floe <Floe>
Output:
        floe with updated fields
"""
function update_floe!(floe)
    collision_force = [0.0 0.0]
    collision_torque = 0.0
    area = floe.area
    mass = floe.mass
    height = floe.height
    intertial_moment = floe.moment

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
                              sin(α)*p[1] + cos(α)p[2]], floe.coords)]

    # Update ice velocities with forces and torques
    dudt = (floe.FxOA*area + collision_force[1])/mass
    dvdt = (floe.FyOA*area + collision_force[2])/mass
    
    frac = if abs(Δt*dudt) > (height/2) && abs(Δt*dvdt) > (height/2)
        frac1 = (sign(dudt)*height/2Δt)/dudt
        frac2 = (sign(dvdt)*height/2Δt)/dvdt
        min(frac1, frac2)
    elseif abs(Δt*dudt) > (height/2) && abs(Δt*dvdt) < (height/2)
        (sign(dudt)*height/2Δt)/dudt
    elseif abs(Δt*dudt) < (height/2) && abs(Δt*dvdt) > (height/2)
        (sign(dvdt)*height/2Δt)/dvdt
    else
        1
    end
    dudt = frac*dudt
    dvdt = frac*dvdt

    floe.u += 1.5Δt*dudt-0.5Δt*floe.p_dudt
    floe.v += 1.5Δt*dvdt-0.5Δt*floe.p_dvdt
    floe.p_dudt = dudt
    floe.p_dvdt = dvdt

    dξdt = (floe.torqueOA*area+collision_torque)/intertial_moment
    dξdt = frac*dξdt
    ξ = floe.ξ + 1.5Δt*dξdt-0.5Δt*floe.p_dξdt
    if abs(ξ) > 1e-5
        ξ = sign(ξ) * 1e-5
    end
    floe.ξ = ξ
    floe.p_dξdt = dξdt

    # TODO: Thermodynamic growth 
    # Questions: What to do about heat flux being a matrix vs float
    # can we just update this at the end of this function, or does it
    # need to happen directly before OA calculations?
    # Calc_trajectory lines 68-73

    # TODO: Floe strain - Calc_trajectory lines 216-288

    # do I need to return floe? Distributed vs not?
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

function run!(simulation)
    println("Model running!")
end