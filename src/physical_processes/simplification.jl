"""
Functions to simplify and reduce individual floes and floe list. 
"""

"""
    dissolve_floe(floe, grid, dissolved)

Dissolve given floe into dissolved ocean matrix.
Inputs:
    floe        <Union{Floe, LazyRow{Floe}}> single floe
    grid        <RegRectilinearGrid> model's grid
    domain      <Domain> model's domain
    dissolve    <Matrix{AbstractFloat}> ocean's dissolved field
Outputs:
    None. Update dissolve matrix with given floe's mass and mark floe for
    removal.
"""
function dissolve_floe!(floe, grid::RegRectilinearGrid, domain, dissolve)
    xidx, yidx = find_grid_cell_index(
        floe.centroid[1],
        floe.centroid[2],
        grid,
    )
    # Wrap indexes around grid if bounds are periodic 
    xidx = shift_cell_idx(xidx, grid.dims[2] + 1, domain.east)
    yidx = shift_cell_idx(yidx, grid.dims[1] + 1, domain.north)
    # If centroid is within bounds after wrapping, add mass to dissolve
    if 0 < xidx <= grid.dims[2] && 0 < yidx <= grid.dims[1]
        dissolve[yidx, xidx] += floe.mass
    end
    floe.status.tag = remove
    return
end

"""
    replace_floe!(
        floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
        new_poly,
        new_mass,
        consts,
        mc_n,
        rng,
    )
Updates existing floe shape and related physical properties based of the polygon
defining the floe.
Inputs:
    floe        <Union{Floe, LazyRow{Floe}}> floe to update
    new_poly    <LG.Polygon> polygon representing new outline of floe
    new_mass    <AbstractFloat> mass of floe
    consts      <Constants> simulation's constants
    mc_n        <Int> number of monte carlo points to attempt to generate
    rng         <RNG> random number generator
Ouputs:
    Updates a given floe's physical properties given new shape and total mass.
"""
function replace_floe!(
    floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    new_poly,
    new_mass,
    consts,
    mc_n,
    rng,
) where {FT}
    # Floe shape
    floe.centroid = find_poly_centroid(new_poly)
    floe.coords = rmholes(find_poly_coords(new_poly)::PolyVec{FT})
    floe.area = LG.area(new_poly)
    floe.height = new_mass/(floe.area * consts.ρi)
    floe.mass = new_mass
    floe.moment = calc_moment_inertia(
        floe.coords,
        floe.centroid,
        floe.height;
        ρi = consts.ρi,
    )
    floe.angles = calc_poly_angles(floe.coords)
    floe.α = FT(0)
    translate!(floe.coords, -floe.centroid[1], -floe.centroid[2])
    floe.rmax = sqrt(maximum([sum(c.^2) for c in floe.coords[1]]))
    # Floe monte carlo points
    mc_x, mc_y, status = generate_mc_points(
        mc_n,
        floe.coords,
        floe.rmax,
        floe.area,
        floe.status,
        rng,
    )
    translate!(floe.coords, floe.centroid[1], floe.centroid[2])
    floe.mc_x = mc_x
    floe.mc_y = mc_y
    # Floe status / identification
    floe.status = status
end

"""
    smooth_floes!(
        floes,
        topography,
        simp_settings,
        collision_settings,
        coupling_settings,
        consts,
        rng,
    )
Smooths floe coordinates for floes with more vertices than the maximum
allowed number. Uses Ramer–Douglas–Peucker algorithm with a user-defined
tolerance. If new shape causes overlap greater with another floe greater than
the maximum percentage allowed, mark the two floes for fusion.
Inputs:
    floes               <StructArray{Floe}> model's floes
    topography          <StructArray{TopographyElement}> domain's topography
    simp_settings       <SimplificationSettings> simulation's simplification
                            settings
    collision_settings  <CollisionSettings> simulation's collision settings
    coupling_settings   <CouplingSettings> simulation's coupling settings
    consts              <Constants> simulation's constants
    rng                 <RNG> random number generator for new monte carlo points
"""
function smooth_floes!(
    floes::StructArray{Floe{FT}},
    topography,
    simp_settings,
    collision_settings,
    coupling_settings,
    consts,
    rng,
) where {FT <: AbstractFloat}
    topo_coords = topography.coords
    Threads.@threads for i in eachindex(floes)
        if length(floes.coords[i][1]) > simp_settings.max_vertices
            poly = LG.simplify(LG.Polygon(floes.coords[i]), simp_settings.tol)
            if !isempty(topo_coords)
                poly = LG.difference(simp_poly, LG.MultiPolygon(topo_coords))
            end
            poly_list = get_polygons(rmholes(poly))::Vector{LG.Polygon}
            new_poly =
                if length(poly_list) == 1
                    poly_list[1]
                else
                    areas = [LG.area(p) for p in poly_list]
                    _, max_idx = findmax(areas)
                    poly_list[max_idx]
                end
            replace_floe!(
                LazyRow(floes, i),
                new_poly,
                floes.mass[i],
                consts,
                coupling_settings.mc_n,
                rng,
            ) 
            # Mark interactions for fusion
            for j in eachindex(floes)
                if i != j && floes.status[j] != remove && potential_interaction(
                    floes.centroid[i],
                    floes.centroid[j],
                    floes.rmax[i],
                    floes.rmax[j],
                )
                    jpoly = LG.Polygon(floes.coords[j])
                    intersect_area = LG.area(LG.intersection(new_poly, jpoly))
                    if intersect_area/LG.area(jpoly) > collision_settings.floe_floe_max_overlap
                        floes.status[i].tag = fuse
                        push!(floes.status[i].fuse_idx, j)
                    end
                end
            end
        end
    end
    return
end

"""
    fuse_floes!(
        floe1,
        floe2,
        consts,
        Δt,
        coupling_settings,
        max_floe_id,
        rng,
    )
Fuses two floes together if they intersect and replaces the larger of the two
floes with their union. Mass and momentum are conserved.
Inputs:
    floe1               <Union{Floe, LazyRow{Floe}}> first floe
    floe2               <Union{Floe, LazyRow{Floe}}> second floe
    consts              <Constants> simulation's constants
    coupling_settings   <CouplingSettings> simulation's coupling settings
    max_floe_id         <Int> maximum floe ID used yet in simulation
    rng                 <RNG> random number generator
Outputs:
    If floes are not intersecting, no changes. If intersecing, the fused floe
    replaces the larger of the two floes and the smaller floe is marked for
    removal. 
"""
function fuse_floes!(
    floe1,
    floe2,
    consts,
    Δt,
    coupling_settings,
    max_floe_id,
    rng,
)
    # Create new polygon if they fuse
    poly1 = LG.Polygon(floe1.coords)::LG.Polygon
    poly2 = LG.Polygon(floe2.coords)::LG.Polygon
    new_poly_list = get_polygons(LG.union(poly1, poly2))::Vector{LG.Polygon}
    if length(new_poly_list) == 1  # if they fused, they will make one polygon
        new_poly = rmholes(new_poly_list[1])
        mass1 = floe1.mass
        mass2 = floe2.mass
        total_mass = mass1 + mass2
        moment1 = floe1.moment
        moment2 = floe2.moment
        x1, y1 = floe1.centroid
        x2, y2 = floe2.centroid
        new_floe =
            if floe1.area >= floe2.area
                floe2.status.tag = remove
                floe1
            else
                floe1.status.tag = remove
                floe2
            end
        replace_floe!(
            new_floe,
            new_poly,
            total_mass,
            consts,
            coupling_settings.mc_n,
            rng,
        )

        # Conservation of momentum in current timestep values
        new_floe.ξ = (1/new_floe.moment) * (
            # averaged spin momentum
            floe1.ξ * moment1 + floe2.ξ * moment2 +
            # change in orbital momentum balanced by change in ξ
            mass1 * (
                floe1.v * (x1 - new_floe.centroid[1]) +
                floe1.u * (y1 - new_floe.centroid[2])
            ) +
            mass2 * (
                floe2.v * (x2 - new_floe.centroid[1]) +
                floe2.u * (y2 - new_floe.centroid[2])
            )
        )
        new_floe.u = (floe1.u * mass1 + floe2.u * mass2)/total_mass
        new_floe.v = (floe1.v * mass1 + floe2.v * mass2)/total_mass

        #= Conservation of momentum in previous timestep values -
            assumes that ratios of floe masses and moments of intertia do not
            change between current and previous timestep. Also assumes simplifed
            timestep equations (i.e. u = p_dxdt + Δt * p_dudt rather than
            u = p_dxdt + 1.5p_dudt - 0.5pp_dudt where pp_dudt is the previous
            value of p_dudt which we no longer have availible as it is replaced
            to give p_dudt)
        =#
        new_p_dxdt = (floe1.p_dxdt * mass1 + floe2.p_dxdt * mass2)/total_mass
        new_p_dydt = (floe1.p_dydt * mass1 + floe2.p_dydt * mass2)/total_mass
        # If new floe existed in previous timestep, estimated centroid location
        new_p_x = new_floe.centroid[1] - Δt * new_p_dxdt
        new_p_y = new_floe.centroid[2] - Δt * new_p_dydt
        new_floe.p_dαdt = (1/new_floe.moment) * (
            # averaged spin momentum in previous timestep
            floe1.p_dαdt * moment1 + floe2.p_dαdt * moment2 +
            # change in orbital momentum balanced by change in previous ξ
            mass1 * (
                floe1.p_dydt * (x1 - Δt * floe1.p_dxdt - new_p_x) +
                floe1.p_dxdt * (y1 - Δt * floe1.p_dydt - new_p_y)
            ) +
            mass2 * (
                floe2.p_dydt * (x2 - Δt * floe2.p_dxdt - new_p_x) +
                floe2.p_dxdt * (y2 - Δt * floe2.p_dydt - new_p_y)
            )
        )   
        new_floe.p_dxdt = new_p_dxdt
        new_floe.p_dydt = new_p_dydt
        # Changes in acceleration in previous timestep determined by velocities
        new_floe.p_dξdt = (new_floe.ξ - new_floe.p_dαdt) / Δt
        new_floe.p_dudt = (new_floe.u - new_floe.p_dxdt) / Δt
        new_floe.p_dvdt = (new_floe.v - new_floe.p_dydt) / Δt

        # Update stress history
        new_floe.stress .= (1/total_mass) * (
            floe1.stress * mass1 .+
            floe2.stress * mass2
        )
        new_floe.stress_history.cb .= (1/total_mass) * (
            floe1.stress_history.cb * mass1 .+
            floe2.stress_history.cb * mass2
        )
        new_floe.stress_history.total = (1/total_mass) * (
            floe1.stress_history.total * mass1 +
            floe2.stress_history.total * mass2)
        # Update IDs
        push!(new_floe.parent_ids, floe1.id)
        push!(new_floe.parent_ids, floe2.id)
        max_floe_id += 1
        new_floe.id = max_floe_id
    end
    return max_floe_id
end

"""
    simplify_floes!(
        model,
        simp_settings,
        collision_settings,
        coupling_settings,
        consts,
        rng,
    )
Simplify the floe list be smoothing vertices, fusing floes, dissolving floes,
and removing floes as needed. 
Inputs:
    model               <Model> model
    simp_settings       <SimplificationSettings> simulation's simplification
                            settings
    collision_settings  <CollisionSettings> simulation's collision settings
    coupling_settings   <CouplingSettings>  simulation's coupling settings
    consts              <Constants> simulation's constants
    rng                 <RNG> random number generator
Outputs:
    Updates floe list and removes floe that won't continue to the next timestep
"""
function simplify_floes!(
    floes,
    model,
    simp_settings,
    collision_settings,
    coupling_settings,
    Δt,
    consts,
    rng,
)
    # Smooth coordinates to reduce total number of vertices
    if smooth_vertices_on
        smooth_floes!(
            floes,
            model.domain.topography,
            simp_settings,
            collision_settings,
            coupling_settings,
            consts,
            rng,
        )
    end
    # Fuse floes that have been marked for fusion
    max_floe_id = model.max_floe_id
    for i in eachindex(floes)
        if floes.status[i] == fuse
            for j in eachindex(floes.status[i].floeidx)
                max_floe_id = fuse_floes!(
                    LazyRow(floes, i),
                    LazyRow(floes, j),
                    consts,
                    Δt,
                    coupling_settings,
                    max_floe_id,
                    rng,
                )
            end
        end
    end
    model.max_floe_id = max_floe_id

    for i in reverse(eachindex(floes))
        if simp_settings.dissolve_on &&
            floes.area[i] < simp_settings.min_floe_area &&
            floes.status[i].tag != remove
            # Dissolve small floes and add mass to ocean
            dissolve_floe!(
                LazyRow(floes, i),
                model.grid,
                model.ocean.dissolved,
            )
            # remove dissolved floes
            StructArrays.foreachfield(
                    field -> deleteat!(field, i),
                    floes,
            )
        # maybe seperate this out??
        elseif floes.status[i].tag == remove
            StructArrays.foreachfield(
                    field -> deleteat!(field, i),
                    floes,
            )
        end
    end
    return
end