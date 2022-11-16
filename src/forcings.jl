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
function calc_OA_forcings!(floe, m)
    c = m.consts
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
    vice = floe.v .+ lx*floe.ξ

    # Force on ice from atmopshere
    uatm = m.wind.u[idx]
    vatm = m.wind.v[idx]
    fx_atm = (c.ρa * c.Cd_ia * sqrt.(uatm.^2 + vatm.^2) .* uatm) .* areas
    fy_atm = (c.ρa * c.Cd_ia * sqrt.(uatm.^2 + vatm.^2) .* vatm) .* areas

    # Force on ice from pressure gradient
    fx_pressure∇ = -ma_ratio * c.f .* m.ocean.v[idx] .* areas
    fy_pressure∇ = ma_ratio * c.f .* m.ocean.u[idx] .* areas

    # Force on ice from ocean
    Δu_OI = m.ocean.u[idx] .- uice
    Δv_OI = m.ocean.v[idx] .- vice
    τx_ocn = c.ρo*c.Cd_io*sqrt.(Δu_OI.^2 + Δv_OI.^2) .* (cos(c.turnθ) .* Δu_OI .- sin(c.turnθ) * Δv_OI)
    τy_ocn = c.ρo*c.Cd_io*sqrt.(Δu_OI.^2 + Δv_OI.^2) .* (sin(c.turnθ) .* Δu_OI .+ cos(c.turnθ) * Δv_OI)
    fx_ocn = τx_ocn .* areas
    fy_ocn = τy_ocn .* areas

    # Sum above forces and find torque
    fx = fx_atm .+ fx_pressure∇ .+ fx_ocn
    fy = fy_atm .+ fy_pressure∇ .+ fy_ocn
    trq = lx.*fy .- ly.*fx  # are these signs right?

    # Add coriolis force to total foces
    fx .+= ma_ratio * c.f * floe.v * areas
    fy .-= ma_ratio * c.f * floe.u * areas

    # Sum forces on ice floe
    floe.fxOA = sum(fx)
    floe.fyOA = sum(fy)
    floe.torqueOA = sum(trq)

    # TODO: Not thread safe
    # Update ocean stress fields with ice on ocean stress
    m.ocean.τx[idx] .= m.ocean.τx[idx].*(1 .- area_ratios) .- τx_ocn.*area_ratios
    m.ocean.τy[idx] .= m.ocean.τy[idx].*(1 .- area_ratios) .- τy_ocn.*area_ratios

    # Update sea-ice fraction
    m.ocean.si_frac[idx] .+= area_ratios
    return
end

function collision_normal_force(c1, c2, region, area, ipoints, force_factor, T)
    force_dir = zeros(T, 2)
    coords = LG.GeoInterface.coordinates(region)::PolyVec{Float64}
    n_ipoints = size(ipoints, 1)
    verts = zeros(Int64, n_ipoints)
    dists = zeros(n_ipoints)
    for i in 1:n_ipoints
        p_i = repeat(ipoints[i, :], length(coords)) # this might be wrong
        dists[i], verts[i] = findmin([sum((c .- p_i).^2) for c in coords[1]])
    end
    dists = sqrt.(dists)
    p = coords[1][verts[findall(d -> d<1, dists)]] # maybe do rmax/1000
    m = length(p)

    # Calculate normal forces
    Δl = zeros(T, 1)
    if m == 2
        Δx = p[2][1] - p[1][1]
        Δy = p[2][2] - p[1][2]
        mag = sqrt(Δx^2 + Δy^2)
        Δl = mag
        if Δl > 0.1  # should match our scale
            force_dir = [-Δy/Δl; Δx/Δl]
        end
    elseif m != 0  # Unusual number of contact points
        x, y = seperate_xy(coords)
        Δx = diff(x)
        xmid = (x[2:end] .+ x[1:end-1]) ./ 2
        Δy = diff(y)
        ymid = (y[2:end] .+ y[1:end-1]) ./ 2
        mag = sqrt.(Δx.^2 .+ Δy.^2)  # vector magniture
        uvec = [-Δy./mag Δx./mag]  # unit vector
        xt = xmid.+uvec[:, 1]./100
        yt = ymid+uvec[:, 1]./100 #  maybe we need to scale?
        in_idx = inpoly2(hcat(xt, yt), hcat(x, y))
        uvec[in_idx[:, 1] .|  in_idx[:, 2], :] *= -1
        Fn = -force_factor * (mag * ones(1, 2)) .* uvec
        dmin_lst = calc_point_poly_dist(xmid, ymid, c1)
        on_idx = findall(d->abs(d)<1e-8, dmin_lst)
        if 0 < length(on_idx) < length(dmin_lst)
            Δl = mean(mag[on_idx])
            if Δl > 0.1
                Fn_tot = sum(Fn[on_idx, :], dims = 1)
                force_dir = Fn_tot'/sqrt(Fn_tot*Fn_tot')
            end
        end
    end
    # Find the direction of the force
    if Δl > 0.1
        c1new = [[c .+ vec(force_dir) for c in c1[1]]]

        # Floe/boudary intersection after being moved in force direction
        new_inter_floe = LG.intersection(LG.Polygon(c1new), LG.Polygon(c2))

        # Need to find which new overlap region corresponds to region k
        new_region_overlaps = LG.getGeometries(LG.intersection(new_inter_floe, region))
        if !isempty(new_region_overlaps)
            ~, max_idx = findmax(LG.area, new_region_overlaps)
            new_overlap_area = LG.area(LG.getGeometries(new_inter_floe)[max_idx])
            if new_overlap_area/area > 1  # Force increased overlap
                force_dir *= -1 
            end
        end
    end
    return (force_dir * area * force_factor)'
end

function calc_collision_forces(c1, c2, regions, region_areas, force_factor, consts, T)
    ipoints = intersect_lines(c1, c2)  # Intersection points
    if isempty(ipoints) || size(ipoints,2) < 2  # No overlap points
        return zeros(T, 1, 2), zeros(T, 1, 2), zeros(T, 1)  # Force, contact points, overlap area
    else
        # Find overlapping regions greater than minumum area
        n1 = length(c1[1]) - 1
        n2 = length(c2[1]) - 1
        min_area = min(n1, n2) * 100 / 1.75
        regions = regions[region_areas .> min_area]
        region_areas = region_areas[region_areas .> min_area]
        overlap = region_areas
        ncontact = length(regions)
        # Calculate forces for each remaining region
        force = zeros(T, ncontact, 2)  # Add type
        pcontact = zeros(T, ncontact, 2)
        for k in eachindex(ncontact)
            normal_force = zeros(T, 1, 2)
            if region_areas[k] == 0
                pcontact[k, :] = zeros(T, 2)
            else
                cx, cy = LG.GeoInterface.coordinates(LG.centroid(regions[k]))::Vector{Float64}
                pcontact[k, :] = [cx, cy]
                normal_force = collision_normal_force(c1, c2, regions[k], region_areas[k], ipoints, force_factor, T)
            end
            # TODO add frictional forces
            G = consts.E/(2*(1+consts.ν))  # Sheer modulus
            friction_force = zeros(T, 1, 2)
            force[k, :] = normal_force .+ friction_force
        end
        return force, pcontact, overlap
    end
end