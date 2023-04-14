"""
Functions to simplify and reduce individual floes and floe list. 
"""

"""
    dissolve_floe(floe, grid, dissolved)

Dissolve given floe into dissolved ocean matrix.
Inputs:
    floe        <Union{Floe, LazyRow{Floe}}>
    grid        <RegRectilinearGrid> model grid
    dissolved   <Matrix{AbstractFloat}> ocean's dissolved field
"""
function dissolve_floe!(floe, grid::RegRectilinearGrid, dissolved)
    xidx, yidx = find_grid_cell_index(
        floe.centroid[1],
        floe.centroid[2],
        grid,
    )
    dissolve[yidx, xidx] += floe.area
    floe.status.tag = remove
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
    coupling_settings,
    max_floe_id,
    rng,
)
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
        centroid1 = floe1.centroid
        centroid2 = floe2.centroid
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
        # Conservation of momentum
        new_floe.u = (floe1.u * mass1 + floe2.u * mass2)/total_mass
        new_floe.v = (floe1.v * mass1 + floe2.v * mass2)/total_mass
        new_floe.ξ = (floe1.ξ * moment1 + floe2.ξ * moment2)/new_floe.moment
        new_floe.p_dudt = (floe1.p_dudt * mass1 + floe2.p_dudt * mass2)/total_mass
        new_floe.p_dvdt = (floe1.p_dvdt * mass1 + floe2.p_dvdt * mass2)/total_mass
        new_floe.p_dxdt = (floe1.p_dxdt * mass1 + floe2.p_dxdt * mass2)/total_mass
        new_floe.p_dydt = (floe1.p_dydt * mass1 + floe2.p_dydt * mass2)/total_mass
        new_floe.p_dξdt = (floe1.p_dξdt * moment1 + floe2.p_dξdt * moment2)/new_floe.moment
        # Update stress history
        new_floe.stress .= (floe1.stress * mass1 .+ floe2.stress * mass2)/total_mass
        new_floe.stress_history.cb .= (floe1.stress_history.cb * mass1 .+
            floe2.stress_history.cb * mass2)/total_mass
        new_floe.stress_history.total = (floe1.stress_history.total * mass1 +
            floe2.stress_history.total * mass2)/total_mass
        # Update IDs
        empty!(new_floe.parent_ids)
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