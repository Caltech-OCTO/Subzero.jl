"""
Functions needed for coupling the ice, ocean, and atmosphere
"""

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
function find_bounding_idx(points, point_idx, gpoints, Δd::Int, end_idx::Int, ::PeriodicBoundary, ::Type{T} = Float64) where T 
    imin, imax = extrema(point_idx)
    Δg = gpoints[end] - gpoints[1]
    knot_idx, knots =
        if imin < 1 && imax > end_idx
            @warn "A floe longer than the domain passed through a periodic boundary. It was removed to prevent overlap."
            Vector{T}(undef, 0)
        else
            imin -= Δd
            imax += (Δd + 1)
            if imin < 1 && imax > end_idx
                [(imin + end_idx):end_idx; 1:end_idx; 1:(imax-end_idx)],
                [gpoints[(imin + end_idx):end_idx] .- Δg; gpoints[1:end_idx]; gpoints[1:(imax-end_idx)] .+ Δg]
            elseif imin < 1
                [(imin + end_idx):end_idx; 1:imax],
                [gpoints[(imin + end_idx):end_idx] .- Δg; gpoints[1:imax]]
            elseif imax > end_idx
                [imin:end_idx; 1:(imax-end_idx)],
                [gpoints[imin:end_idx]; gpoints[1:(imax-end_idx)] .+ Δg]
            else
                [imin:imax;]
            end
        end
    return knots, knot_idx, points, point_idx
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
function find_bounding_idx(points, point_idx, gpoints, Δd::Int, end_idx::Int, ::AbstractBoundary, ::Type{T} = Float64) where T
    points = points[point_idx .> 0]
    point_idx = point_idx[point_idx .> 0]
    imin, imax = extrema(point_idx)
    imin -= Δd
    imax += (Δd + 1) # point in ith gridcell is between the i and i+1 grid line
    imin = (imin < 1) ? 1 : imin
    imax = (imax > end_idx) ? end_idx : imax
    return gpoints[imin:imax], [imin:imax;], points, point_idx
end

function point_grid_indices(xp, yp, grid::Grid)
    xidx = floor.(Int, (xp .- grid.xg[1])/(grid.xg[2] - grid.xg[1])) .+ 1
    yidx = floor.(Int, (yp .- grid.yg[1])/(grid.yg[2] - grid.yg[1])) .+ 1
    return xidx, yidx
end

"""
    floe_OA_forcings!(floe, m, t::Type{T} = Float64)

Calculate the effects on the ocean and atmpshere on floe i within the given model
and the effects of the ice floe on the ocean.

Inputs:
        floe    <Floe> floe
        m       <Model> given model
        c       <Constants> constants within simulation
        Δd      <Int> buffer for monte carlo interpolation knots
Outputs:
        None. Both floe and ocean fields are updated in-place.
"""
function floe_OA_forcings!(floe, m, c, Δd, ::Type{T} = Float64) where T
    # Rotate and translate Monte Carlo points to current floe location
    α = floe.α
    α_rot = [cos(α) -sin(α); sin(α) cos(α)]
    mc_p = α_rot * [floe.mc_x'; floe.mc_y']
    mc_xr = mc_p[1, :] .+ floe.centroid[1]
    mc_yr = mc_p[2, :] .+ floe.centroid[2]
    mc_xidx, mc_yidx = point_grid_indices(mc_xr, mc_yr, m.grid)

    nrows, ncols = m.grid.dims
    xknots, xknot_idx, mc_xr, mc_xidx = find_bounding_idx(mc_xr, mc_xidx, m.grid.xg, Δd, ncols, m.domain.east, T)
    yknots, yknot_idx, mc_yr, mc_yidx = find_bounding_idx(mc_yr, mc_yidx, m.grid.yg, Δd, nrows, m.domain.north, T)

    if isempty(xknots) && isempty(yknots)
        floe.alive = 0
    else
        # Wind Interpolation for Monte Carlo Points
        uatm_interp = linear_interpolation((xknots, yknots), m.wind.u[xknot_idx, yknot_idx])
        vatm_interp = linear_interpolation((xknots, yknots), m.wind.v[xknot_idx, yknot_idx])
        avg_uatm = mean([uatm_interp(mc_xr[i], mc_yr[i]) for i in eachindex(mc_xr)])
        avg_vatm = mean([vatm_interp(mc_xr[i], mc_yr[i]) for i in eachindex(mc_xr)])

        # Ocean Interpolation for Monte Carlo Points
        uocn_interp = linear_interpolation((xknots, yknots), m.ocean.u[xknot_idx, yknot_idx])
        vocn_interp = linear_interpolation((xknots, yknots), m.ocean.v[xknot_idx, yknot_idx])
        hflx_interp = linear_interpolation((xknots, yknots), m.ocean.v[xknot_idx, yknot_idx])
        uocn = [uocn_interp(mc_xr[i], mc_yr[i]) for i in eachindex(mc_xr)]
        vocn = [vocn_interp(mc_xr[i], mc_yr[i]) for i in eachindex(mc_xr)]
        floe.hflx = mean([hflx_interp(mc_xr[i], mc_yr[i]) for i in eachindex(mc_xr)])

        # Force on ice from atmopshere
        τx_atm = (c.ρa * c.Cd_ia * sqrt(avg_uatm^2 + avg_vatm^2) * avg_uatm)
        τy_atm = (c.ρa * c.Cd_ia * sqrt(avg_uatm^2 + avg_vatm^2) * avg_vatm)

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