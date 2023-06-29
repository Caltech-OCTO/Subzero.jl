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
    dissolved    <Matrix{AbstractFloat}> ocean's dissolved field
Outputs:
    None. Update dissolved matrix with given floe's mass and mark floe for
    removal.
"""
function dissolve_floe!(floe, grid::RegRectilinearGrid, domain, dissolved)
    xidx, yidx = find_grid_cell_index(
        floe.centroid[1],
        floe.centroid[2],
        grid,
    )
    # Wrap indexes around grid if bounds are periodic 
    xidx = shift_cell_idx(xidx, grid.Nx + 1, domain.east)
    yidx = shift_cell_idx(yidx, grid.Ny + 1, domain.north)
    # If centroid is within bounds after wrapping, add mass to dissolved
    if 0 < xidx <= grid.Nx && 0 < yidx <= grid.Ny
        dissolved[yidx, xidx] += floe.mass
    end
    return
end

"""
    smooth_floes!(
        floes,
        topography,
        simp_settings,
        collision_settings,
        coupling_settings,
        consts,
        Δt,
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
    Δt                  <Int> length of simulation timestep in seconds
    rng                 <RNG> random number generator for new monte carlo points
"""
function smooth_floes!(
    floes::StructArray{Floe{FT}},
    topography,
    simp_settings,
    collision_settings,
    coupling_settings,
    consts,
    Δt,
    rng,
) where {FT <: AbstractFloat}
    topo_coords = topography.coords
    for i in eachindex(floes)
        if length(floes.coords[i][1]) > simp_settings.max_vertices
            poly = LG.simplify(LG.Polygon(floes.coords[i]), simp_settings.tol)
            if !isempty(topo_coords)
                poly = LG.difference(poly, LG.MultiPolygon(topo_coords))
            end
            poly_list = get_polygons(rmholes(poly))::Vector{LG.Polygon}
            simp_poly =
                if length(poly_list) == 1
                    poly_list[1]
                else
                    areas = [LG.area(p) for p in poly_list]
                    _, max_idx = findmax(areas)
                    poly_list[max_idx]
                end
            x_tmp, y_tmp = floes.centroid[i]
            moment_tmp = floes.moment[i]
            replace_floe!(
                LazyRow(floes, i),
                simp_poly,
                floes.mass[i],
                consts,
                coupling_settings.npoints,
                rng,
            )
            # conserve momentum
            conserve_momentum_combination!(
                floes.mass[i],
                moment_tmp,
                x_tmp,
                y_tmp,
                Δt,
                LazyRow(floes, i),
            )
            # Mark interactions for fusion
            for j in eachindex(floes)
                if i != j && floes.status[j].tag != remove && potential_interaction(
                    floes.centroid[i],
                    floes.centroid[j],
                    floes.rmax[i],
                    floes.rmax[j],
                )
                    if floes.status[j].tag == fuse && i in floes.status[j].fuse_idx
                        floes.status[i].tag = fuse
                        push!(floes.status[i].fuse_idx, j)
                    else
                        jpoly = LG.Polygon(floes.coords[j])
                        intersect_area = LG.area(LG.intersection(simp_poly, jpoly))
                        if intersect_area/LG.area(jpoly) > collision_settings.floe_floe_max_overlap
                            floes.status[i].tag = fuse
                            push!(floes.status[i].fuse_idx, j)
                        end
                    end
                end
            end
        end
    end
    return
end

"""
    fuse_two_floes!(
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
function fuse_two_floes!(
    floe1,
    floe2,
    consts,
    Δt,
    coupling_settings,
    max_floe_id,
    prefuse_max_floe_id,
    rng,
)
    # Create new polygon if they fuse
    poly1 = LG.Polygon(floe1.coords)::LG.Polygon
    poly2 = LG.Polygon(floe2.coords)::LG.Polygon
    new_poly_list = get_polygons(LG.union(poly1, poly2))::Vector{LG.Polygon}
    if length(new_poly_list) == 1  # if they fused, they will make one polygon
        new_poly = rmholes(new_poly_list[1])
        keep_floe, remove_floe =
            if floe1.area >= floe2.area
                floe1, floe2
            else
                floe2, floe1
            end
        remove_floe.status.tag = remove
        # record as value will change with replace
        mass_tmp = keep_floe.mass
        moment_tmp = keep_floe.moment
        x_tmp, y_tmp = keep_floe.centroid
        # give keep_floe new shape
        replace_floe!(
            keep_floe,
            new_poly,
            floe1.mass + floe2.mass,
            consts,
            coupling_settings.npoints,
            rng,
        )
        # conserve momentum
        conserve_momentum_combination!(
            mass_tmp,
            moment_tmp,
            x_tmp,
            y_tmp,
            Δt,
            keep_floe,
            remove_floe,
        )
        # Update stress history
        keep_floe.stress .= (1/keep_floe.mass) * (
            keep_floe.stress * mass_tmp .+
            remove_floe.stress * remove_floe.mass
        )
        keep_floe.stress_history.cb .= (1/keep_floe.mass) * (
            keep_floe.stress_history.cb * mass_tmp .+
            remove_floe.stress_history.cb * remove_floe.mass
        )
        keep_floe.stress_history.total = (1/keep_floe.mass) * (
            keep_floe.stress_history.total * mass_tmp +
            remove_floe.stress_history.total * remove_floe.mass)
        # Update IDs
        if floe1.id <= prefuse_max_floe_id
            push!(keep_floe.parent_ids, floe1.id)
        end
        if floe2.id <= prefuse_max_floe_id
            push!(keep_floe.parent_ids, floe2.id)
        end
        max_floe_id += 1
        keep_floe.id = max_floe_id
    end
    return max_floe_id
end

"""
    fuse_floes!(
        floes,
        max_floe_id,
        coupling_settings,
        Δt,
        consts,
        rng,
    )

Fuse all floes marked for fusion.
Inputs:
    floes               <StructArray{Floe}> model's floes
    max_floe_id         <Int> maximum floe ID created yet
    coupling_settings   <CouplingSettings>  simulation's coupling settings
    Δt                  <Int> simulation timestep in seconds
    consts              <Constants> simulation's constants
    rng                 <RNG> random number generator
Outputs:
    None. Fuses floes marked for fusion. Marks floes fused into another floe
    for removal. 
"""
function fuse_floes!(
    floes,
    max_floe_id,
    coupling_settings,
    Δt,
    consts,
    rng,
)
    prefuse_max_floe_id = max_floe_id
    for i in eachindex(floes)
        if floes.status[i].tag == fuse
            for j in floes.status[i].fuse_idx
                # not already fused with another floe or marked for removal
                if floes.status[j].tag != remove
                    max_floe_id = fuse_two_floes!(
                        LazyRow(floes, i),
                        LazyRow(floes, j),
                        consts,
                        Δt,
                        coupling_settings,
                        max_floe_id,
                        prefuse_max_floe_id,
                        rng,
                    )
                end
            end
        end
    end
    return max_floe_id
end

"""
    remove_floes!(
        floes,
        grid,
        domain,
        dissolved,
        simp_settings
    )

Remove floes marked for removal and dissolve floes smaller than minimum floe
area if the dissolve setting is on.
Inputs:
    floes           <StructArray{Floe}> model's floes
    grid            <AbstractGrid> model's grid
    domain          <Domain> model's domain
    dissolved       <Matrix{AbstractFloat}> ocean's dissolved field
    simp_settings   <SimplificationSettings> simplification settings
Outputs:
    None. Removes floes that do not continue to the next timestep and reset all
    continuing floes status to active.
"""
function remove_floes!(
    floes,
    grid,
    domain,
    dissolved,
    simp_settings
)
    for i in reverse(eachindex(floes))
        if floes.status[i].tag != remove && (
            floes.area[i] < simp_settings.min_floe_area ||
            floes.height[i] < simp_settings.min_floe_height
        )
            # Dissolve small/thin floes and add mass to ocean
            dissolve_floe!(
                LazyRow(floes, i),
                grid,
                domain,
                dissolved,
            )
            # remove dissolved floes
            StructArrays.foreachfield(
                    field -> deleteat!(field, i),
                    floes,
            )
        elseif floes.status[i].tag == remove
            StructArrays.foreachfield(
                    field -> deleteat!(field, i),
                    floes,
            )
        else  # reset continuing floes' status
            floes.status[i].tag = active
            empty!(floes.status[i].fuse_idx)
        end
    end
    return
end

"""
    simplify_floes!(
        model,
        simp_settings,
        collision_settings,
        coupling_settings,
        Δt,
        consts,
        rng,
    )
Simplify the floe list be smoothing vertices, fusing floes, dissolving floes,
and removing floes as needed. 
Inputs:
    model               <Model> model
    max_floe_id         <Int> maximum floe id in simulation
    simp_settings       <SimplificationSettings> simulation's simplification
                            settings
    collision_settings  <CollisionSettings> simulation's collision settings
    coupling_settings   <CouplingSettings>  simulation's coupling settings
    Δt                  <Int> simulation timestep in seconds
    consts              <Constants> simulation's constants
    rng                 <RNG> random number generator
Outputs:
    Updates floe list and removes floe that won't continue to the next timestep
"""
function simplify_floes!(
    model,
    max_floe_id,
    simp_settings,
    collision_settings,
    coupling_settings,
    Δt,
    consts,
    rng,
)
    # Smooth coordinates to reduce total number of vertices
    if simp_settings.smooth_vertices_on
        smooth_floes!(
            model.floes,
            model.domain.topography,
            simp_settings,
            collision_settings,
            coupling_settings,
            consts,
            Δt,
            rng,
        )
    end
    # Fuse floes that have been marked for fusion
    max_floe_id = fuse_floes!(
        model.floes,
        max_floe_id,
        coupling_settings,
        Δt,
        consts,
        rng,
    )

    # Remove floes marked for removal and dissolving
    remove_floes!(
        model.floes,
        model.grid,
        model.domain,
        model.ocean.dissolved,
        simp_settings
    )
    return max_floe_id
end