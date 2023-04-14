function find_cell_indices(xp, yp, grid::RegRectilinearGrid)
    xidx = floor.(Int, (xp .- grid.xg[1])/(grid.xg[2] - grid.xg[1])) .+ 1
    yidx = floor.(Int, (yp .- grid.yg[1])/(grid.yg[2] - grid.yg[1])) .+ 1
    return xidx, yidx
end

function dissolve_small_floes(model, simp_settings)
    small_floe_idx = (model.floes.area .< simp_settings.min_floe_area) .&&
        model.floes.active  # if floe is active and small then 'melt' it
    small_floe_centroids = model.floes.centroid[small_floe_idx]
    small_floe_areas = model.floes.area[small_floe_idx]
    xidx, yidx = find_cell_indices(
        first.(small_floe_centroids),
        last.(small_floe_centroids),
        model.grid,
    )
    for i in eachindex(xidx)
        model.ocean.dissolve[yidx[i], xidx[i]] += small_floe_areas[i]
    end
end

function simplify_floes!(
    floes::StructArray{Floe{FT}},
    topography,
    simp_settings,
    collision_settings,
    coupling_settings,
    consts,
    rng,
) where {FT <: AbstractFloat}
    topo_coords = topography.coords
    for i in eachindex(floes)
        poly = LG.simplify(LG.Polygon(floes.coords[i]), simp_settings.tol)
        if !isempty(topo_coords)
            poly = LG.difference(simp_poly, LG.MultiPolygon(topo_coords))
        end
        poly_lst = get_polygons(rmholes(poly))::Vector{LG.Polygon}
        new_poly =
            if length(poly_list) == 1
                poly_list[1]
            else
                areas = [LG.area(p) for p in poly_lst]
                _, max_idx = findmax(areas)
                poly_lst[max_idx]
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
                floes.coords[i],
                floes.coords[j],
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
    return
end

"""
    fuse_floes(
        floes::Vector{Union{Floe{FT}}, LazyRow{Floe{FT}}},
    )

Fuse two floes together into their union. If there is no union, leave original
floes unchanged. 
Inputs:
    floes       <Vector{Floe}> vector of floes to fuse together
Outputs:
    new_floes   <StructArray{Floe}> new fused floe(s) to add to model floe list
"""
# function find_fuse_groups(floes)
#     fuse_groups = Dict{Int, Set{Int}}()
#     for i in eachindex(floes)
#         if floes.status[i] == fuse_floes
#             sort!(floes.status[i].fuse_idx)
#             min_idx = floes.status[i].fuse_idx[1]  # smallest is now first
#             if min_idx > i  # have not seen a fusion group involving this floe
#                 fuse_groups[i] = floes.status[i].fuse_idx
#             else  # this floe is in a fusion group with a lesser-index floe
#                 while !haskey(fuse_groups, min_idx)
#                     min_idx = floes.status[min_idx].fuse_idx[1]
#                 end
#                 # Add new indices to group and delete any duplicate groups
#                 for j in eachindex(floes.status[i].fuse_idx[2:end])
#                     push!(fuse_groups[min_idx], floes.status[i].fuse_idx[j])
#                     delete!(fuse_groups, floes.status[i].fuse_idx[j])
#                 end
#             end
#         end
#     end
# end

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
    poly1 = LG.Polygon(floe1)::LG.Polygon
    poly2 = LG.Polygon(floe2)::LG.Polygon
    new_poly_list = get_polygons(LG.union(poly1, poly2))::Vector{LG.Polygon}
    if length(new_poly_list) == 1  # if they fused, they will make one polygon
        new_poly = rmholes(new_poly_list[1])
        mass1 = floe1.mass
        mass2 = floe2.mass
        moment1 = floe1.moment
        moment2 = floe2.moment
        total_mass = mass1 + mass2
        new_floe =
            if floe1.area > floe2.area
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
        empty!(floe.parent_ids)
        push!(floe.parent_ids, floe1.id)
        push!(floe.parent_ids, floe2.id)
        new_floe.id = max_floe_id
        max_floe_idx += 1
    end
    return max_floe_id
end