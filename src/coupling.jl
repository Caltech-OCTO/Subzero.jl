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
    domain_coords(domain::Domain)
Inputs:
        domain<Domain>
Output:
        RingVec coordinates for edges of rectangular domain based off of boundary values
"""
function cell_coords(xmin, xmax, ymin, ymax)
    return [[[xmin, ymax], [xmin, ymin],
             [xmax, ymin], [xmax, ymax],
             [xmin, ymax]]]
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
    find_bounding_idx(point_idx, Δd::Int, end_idx::Int, ::PeriodicBoundary, ::Type{T} = Float64)

Find indicies in list of grid lines that surround points with indicies 'point_idx'
with a buffer of Δd indices on each side of the points. In this case, the points are being considered near a
periodic boundary, which means that they can loop around to the other side of the grid. 
Inputs:
        point_idx <Vector{Int}> vector of indices representing indices of a list of points on the grid
        Δd        <Int> number of buffer grid cells to include on either side of the provided indicies
        end_idx   <Int> number of grid lines 
Outputs:
        List of indices that include, and surround the given indices with a buffer of Δd on each side.
        Assumes that domain is periodic, so if the buffer causes the indices to be larger or smaller
        than the number of grid lines, the indices will wrap around to the other side of the grid. 
"""
function find_bounding_idx(point_idx, Δd::Int, end_idx::Int, ::PeriodicBoundary, ::Type{T} = Float64) where T
    imin, imax = extrema(point_idx)
    imin -= Δd
    imax += Δd
    bounding =
        if imin < 1 && imax > end_idx
            @warn "A floe longer than the domain passed through a periodic boundary. It was removed to prevent overlap."
            Vector{T}(undef, 0)
        elseif imin < 1
            [1:max_idx; (imin + end_idx):end_idx]
        elseif imax > end_idx
            [1:(imax-end_idx); imin:end_idx]
        else
            [imin:imax;]
        end
    return bounding
end

"""
    find_bounding_idx(point_idx, Δd::Int, end_idx::Int, ::AbstractBoundary, ::Type{T} = Float64)

Find indicies in list of grid lines that surround points with indicies 'point_idx'
with a buffer of Δd indices on each side of the points. In this case, the points are being considered near a
non-periodic boundary, so the points must be within the grid. 
    Inputs:
            point_idx <Vector{Int}> vector of indices representing indices of a list of points on the grid
            Δd        <Int> number of buffer grid cells to include on either side of the provided indicies
            end_idx   <Int> number of grid lines 
    Outputs:
            List of indices that include, and surround the given indices with a buffer of Δd on each side.
            Assumes that domain is non-periodic so the points are cut-off at the edge of the grid,
            even if this means that there is no buffer.
"""
function find_bounding_idx(point_idx, Δd::Int, end_idx::Int, ::AbstractBoundary, ::Type{T} = Float64) where T
    imin, imax = extrema(point_idx)
    imin -= Δd
    imax += Δd
    imin = (imin < 1) ? 1 : imin
    imax = (imax > end_idx) ? end_idx : imax
    return [imin:imax;]
end

"""
    floe_OA_forcings!(floe, m, t::Type{T} = Float64)

Calculate the effects on the ocean and atmpshere on floe i within the given model
and the effects of the ice floe on the ocean.

Inputs:
        floe    <Floe> floe
        m       <Model> given model
Outputs:
        None. Both floe and ocean fields are updated in-place.
Note: For floes that are completly out of the Grid, simulation will error. 
"""
function floe_OA_forcings!(floe, m, c, ::Type{T} = Float64) where T
    # Rotate and translate Monte Carlo points to current floe location
    α = floe.α
    α_rot = [cos(α) -sin(α); sin(α) cos(α)]
    mc_p = α_rot * [floe.mc_x'; floe.mc_y']
    mc_xr = mc_p[1, :] .+ floe.centroid[1]
    mc_yr = mc_p[2, :] .+ floe.centroid[2]

    # Find grid cell that each point is in - indices might not be in bounds if floe is past edge of domain
    mc_xidx = Int.(fld.(mc_xr .- m.grid.xg[1], m.grid.xg[2] - m.grid.xg[1])) .+ 1
    mc_yidx = Int.(fld.(mc_yr .- m.grid.yg[1], m.grid.yg[2] - m.grid.yg[1])) .+ 1

    # Find grid indicies surrounding floe with buffer
    nrows, ncols = m.grid.dims
    xbound_idx = find_bounding_idx(mc_xidx, 2, ncols, m.domain.east, T)
    ybound_idx = find_bounding_idx(mc_yidx, 2, nrows, m.domain.north, T)

    if isempty(xbound_idx) && isempty(ybound_idx)
        floe.alive = 0
    else
        # Adjust indices such that they are all within the domain
        mc_xidx[mc_xidx .< 1] .+= ncols
        mc_xidx[mc_xidx .> ncols] .-= ncols
        mc_yidx[mc_yidx .< 1] .+= nrows
        mc_yidx[mc_yidx .> nrows] .-= nrows

        # Floe heatflux for grid cells below monte carlo points
        floe.hflx = mean(m.ocean.hflx[mc_xidx, mc_yidx])

        # Grid line values for indices surrounding floe
        xg_interp = m.grid.xg[xbound_idx]
        yg_interp = m.grid.yg[ybound_idx]

        # Wind Interpolation for Monte Carlo Points
        uatm_interp = linear_interpolation((xg_interp, yg_interp), m.wind.u[xbound_idx, ybound_idx])
        vatm_interp = linear_interpolation((xg_interp, yg_interp), m.wind.v[xbound_idx, ybound_idx])
        uatm = mean([uatm_interp(mc_xr[i], mc_yr[i]) for i in eachindex(mc_xr)])
        vatm = mean([vatm_interp(mc_xr[i], mc_yr[i]) for i in eachindex(mc_xr)])

        # Ocean Interpolation for Monte Carlo Points
        uocn_interp = linear_interpolation((xg_interp, yg_interp), m.ocean.u[xbound_idx, ybound_idx])
        vocn_interp = linear_interpolation((xg_interp, yg_interp), m.ocean.v[xbound_idx, ybound_idx])
        uocn = [uocn_interp(mc_xr[i], mc_yr[i]) for i in eachindex(mc_xr)]
        vocn = [vocn_interp(mc_xr[i], mc_yr[i]) for i in eachindex(mc_xr)]

        # Force on ice from atmopshere
        τx_atm = (c.ρa * c.Cd_ia * sqrt(uatm^2 + vatm^2) * uatm)
        τy_atm = (c.ρa * c.Cd_ia * sqrt(uatm^2 + vatm^2) * vatm)

        # Force on ice from pressure gradient
        ma_ratio = floe.mass/floe.area
        τx_pressure∇ = -ma_ratio * c.f .* vocn
        τy_pressure∇ = ma_ratio * c.f .* uocn

        # Find total velocity at each monte carlo point
        mc_rad = sqrt.(mc_p[1, :].^2 .+ mc_p[1, :].^2)
        mc_θ = atan.(mc_p[1, :], mc_p[2, :])
        mc_u = floe.u .- floe.ξ * mc_rad .* sin.(mc_θ)
        mc_v = floe.v .+ floe.ξ * mc_rad .* cos.(mc_θ)

        # Force on ice from ocean
        Δu_OI = uocn .- mc_u
        Δv_OI = vocn .- mc_v
        τx_ocn = c.ρo*c.Cd_io*sqrt.(Δu_OI.^2 + Δv_OI.^2) .* (cos(c.turnθ) .* Δu_OI .- sin(c.turnθ) * Δv_OI)
        τy_ocn = c.ρo*c.Cd_io*sqrt.(Δu_OI.^2 + Δv_OI.^2) .* (sin(c.turnθ) .* Δu_OI .+ cos(c.turnθ) * Δv_OI)

        # Save ocean stress fields to update ocean
        τx_group = Dict{Tuple{Int64, Int64}, Vector{T}}()
        τy_group = Dict{Tuple{Int64, Int64}, Vector{T}}()
        for i in eachindex(mc_xidx)
            idx = (mc_xidx[i], mc_yidx[i])
            if haskey(τx_group, idx)
                push!(τx_group[idx], τx_ocn[i])
                push!(τy_group[idx], τy_ocn[i])
            else
                τx_group[idx] = [τx_ocn[i]]
                τy_group[idx] = [τy_ocn[i]]
            end
        end

        floe_poly = LG.Polygon(floe.coords)
        for key in keys(τx_group)
            grid_poly = LG.Polygon(cell_coords(m.grid.xg[key[1]], m.grid.xg[key[1] + 1], m.grid.yg[key[2]], m.grid.yg[key[2] + 1]))
            floe_area_in_cell = LG.area(LG.intersection(grid_poly, floe_poly))
            fx_mean = mean(τx_group[key]) * floe_area_in_cell
            fy_mean = mean(τy_group[key]) * floe_area_in_cell
            # Need to add a lock here...
            m.ocean.τx[key[1], key[2]] += fx_mean
            m.ocean.τy[key[1], key[2]] += fy_mean
            m.ocean.si_area[key[1], key[2]] += floe_area_in_cell
        end


        # Sum above forces and find torque
        τx = τx_atm .+ τx_pressure∇ .+ τx_ocn
        τy = τy_atm .+ τy_pressure∇ .+ τy_ocn
        τtrq = (-τx .* sin.(mc_θ) .+ τy .* cos.(mc_θ)) .* mc_rad

        # Add coriolis force to total foces - does not contribute to torque
        τx .+= ma_ratio * c.f * floe.v
        τy .-= ma_ratio * c.f * floe.u

        # Sum forces on ice floe
        floe.fxOA = mean(τx) * floe.area
        floe.fyOA = mean(τy) * floe.area
        floe.trqOA = mean(τtrq) * floe.area
    end
    return
end