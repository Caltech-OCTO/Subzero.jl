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
    cell_area_ratio(cell_poly, floe_poly)

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
    floe_poly = LG.Polygon(floe.coords)
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
    # y values are rows and x values are columns
    idx = CartesianIndex.(Tuple.(eachrow(hcat(yidx,xidx))))
    return area_ratios, xidx, yidx, idx
end

"""
    calc_OA_forcings!(m, i)

Calculate the effects on the ocean and atmpshere on floe i within the given model
and the effects of the ice floe on the ocean.

Inputs:
        floe    <Floe> floe
        m       <Model> given model
Outputs:
        None. Both floe and ocean fields are updated in-place.
Note: For floes that are completly out of the Grid, simulation will error. 
"""
function floe_OA_forcings!(floe, m)
    c = m.consts
    Δx = m.grid.xg[2] - m.grid.xg[1]
    Δy = m.grid.yg[2] - m.grid.yg[1]

    xi = floe.centroid[1]
    yi = floe.centroid[2]
    α_rot = [cos(floe.α) -sin(α); sin(α) cos(α)]
    # Rotate and translate Monte Carlo points to current floe location
    mc_p = α_rot * [floe.mc_x'; floe.mc_y']
    mc_xr = mc_p[1, :] .+ xi
    mc_xr = mc_p[2, :] .+ yi

    # Find total velocity at each monte carlo point
    mc_rad = sqrt.(mc_p[1, :]^2 .+ mc_p[1, :]^2)
    mc_θ = atan.(mc_p[1, :], mc_p[2, :])
    mc_u = floe.u .- floe.ξ * mc_rad .* sin(mc_θ)
    mc_v = floe.v .+ floe.ξ * mc_rad .* cos(mc_θ)
    mc_xidx = cld(mc_xr .- m.grid.xg[1], Δx)
    mc_yidx = cld(mc_yr .- m.grid.yg[1], Δy)

    # Find grid indicies surrounding floe --> can we use the mc indicies for this??
    rbound = sqrt(2) * floe.rmax
    xbound = rbound + 2Δx
    ybound = rbound + 2Δy
    xbound_idx = -xbound .<= m.grid.xg .- xi .<= xbound
    ybound_idx = -ybound .<= m.grid.yg .- yi .<= ybound
    xg_idx = m.grid.xg[xbound_idx]
    yg_idx = m.grid.xg[ybound_idx]
    # TODO: Interp ocean and wind

    # Floe heatflux
    #floe.hflx = mean(m.hflx[idx])

    # Force on ice from atmopshere
    # TODO: Replace with Interp mean
    uatm = m.wind.u[idx]
    vatm = m.wind.v[idx]
    τx_atm = (c.ρa * c.Cd_ia * sqrt(uatm^2 + vatm^2) * uatm)
    τy_atm = (c.ρa * c.Cd_ia * sqrt(uatm^2 + vatm^2) * vatm)

    # Force on ice from pressure gradient
    # TODO: Replace with Interp 
    uocn = m.ocean.u[idx]
    vocn = m.ocean.v[idx]
    τx_pressure∇ = -ma_ratio * c.f .* vocn
    τy_pressure∇ = ma_ratio * c.f .* uocn

    # Force on ice from ocean
    Δu_OI = uocn .- mc_u
    Δv_OI = vocn .- mc_v
    τx_ocn = c.ρo*c.Cd_io*sqrt.(Δu_OI.^2 + Δv_OI.^2) .* (cos(c.turnθ) .* Δu_OI .- sin(c.turnθ) * Δv_OI)
    τy_ocn = c.ρo*c.Cd_io*sqrt.(Δu_OI.^2 + Δv_OI.^2) .* (sin(c.turnθ) .* Δu_OI .+ cos(c.turnθ) * Δv_OI)

    # Save ocean stress fields to update ocean
    floe.mc_xocnτ[mc_xidx, mc_yidx] .= τx_ocn
    floe.mc_yocnτ[mc_xidx, mc_yidx] .= τy_ocn

    # Sum above forces and find torque
    τx = τx_atm .+ τx_pressure∇ .+ τx_ocn
    τy = τy_atm .+ τy_pressure∇ .+ τy_ocn
    τtrq = (-τx .* sin.(mc_θ) .+ τy .*cos(mc_θ)) .* mc_rad

    # Add coriolis force to total foces
    τx .+= ma_ratio * c.f * floe.v
    τy .-= ma_ratio * c.f * floe.u

    # Sum forces on ice floe
    floe.fxOA = mean(τx) * floe.area
    floe.fyOA = mean(τy) * floe.area
    floe.trqOA = mean(τtrq) * floe.area
    return
end

function update_ocean_stress(ocean)

end