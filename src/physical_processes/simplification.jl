function find_cell_indices(xp, yp, grid::RegRectilinearGrid)
    xidx = floor.(Int, (xp .- grid.xg[1])/(grid.xg[2] - grid.xg[1])) .+ 1
    yidx = floor.(Int, (yp .- grid.yg[1])/(grid.yg[2] - grid.yg[1])) .+ 1
    return xidx, yidx
end

function dissolve_small_floes(model, simp_settings)
    small_floe_idx = (model.floes.area .< simp_settings.min_floe_area) .&&
        model.floes.alive  # if floe is alive and small then 'melt' it
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

function remove_transfer_floes(remove, transfer, n_init_floes, ocean.dissolved)
    transfer = transfer[transfer .<= n_init_floes]
end

function simplify_floe(
    i,
    floes::StructArray{Floe{FT}},
    topo_poly,
    simp_settings,
    collision_settings,
) where {FT <: AbstractFloat}
    poly = LG.simplify(LG.Polygon(floes.coords[i]), simp_settings.tol)
    # Remove overlaps with topography
    if !isempty(domain.topography)
        poly = LG.difference(simp_poly, topo_poly)
    end
    simp_poly = rmholes(poly)
    simp_coords = find_poly_coords(simp_poly)::PolyVec{FT}
    simp_centroid = find_poly_centroid(simp_poly)::Vector{FT}
    # Scale floe by area change
    simp_area = LG.area(simp_poly)::FT
    translate!(simp_coords, [-simp_centroid[1], -simp_centroid[2]])
    simp_coords .*= sqrt(floes.area[i]/simp_area)
    translate!(simp_coords, [simp_centroid[1], simp_centroid[2]])
    simp_poly = LG.Polygon(simp_coords)  # what if this is multiple pieces... --> we will just take the largest??
    for j in eachindex(floes)
        if i != j && (sum((floes.centroid[i] .- floes.centroid[j]).^2) < (floes.rmax[i] + floes.rmax[j])^2)
            intersect_poly = LG.Polygon(floes.coords[j])
            intersect_area = LG.area(LG.intersection(simp_poly, intersect_poly))
            if intersect_area/floes.area[j] > collision_settings.floe_floe_max_overlap
                # set status to fuse
            end
        end
    end
    # re-do all of the floe's stuff


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

#=
Let's use the new idea on the board. Let's take in all floes and go through
and create the groupings.

Issues:
For collisions, the floes have since moved so they might not be overlapping and
the union could return two seperate floes
For simplification, this won't happen since we know that floe's are still in the
same spot...
=#
function fuse_floes(
    floes::Vector{Union{Floe{FT}}, LazyRow{Floe{FT}}},
) where {FT <: AbstractFloat}
    polys = [LG.Polygon(f) for f in floes]
    new_poly = reduce(LG.union, polys)
    new_poly_list = LG.getGeometries(new_poly)
    for i in eachindex(new_poly_list)
        new_poly_list[i] = rmholes(new_poly_list[i])
    end
    new_floes = StructArray{Floe{FT}}(undef, 0)
    total_mass = FT(0)
    id_lst = Vector{Int}()
    for f in floes
        total_mass += f.mass
        push!(id_lst, f.id)
    end
    create_floes_conserve_mass!(  # Do we want to just assume one floe??
        new_poly_list,
        total_mass::FT,
        u::FT, 
        v, 
        ξ,
        mc_n,
        nhistory,
        rng,
        consts,
        new_floes,
    )
    if !isempty(new_floes)
        centroid = find_poly_centroid(new_poly)::Vector{FT}
        total_inertia = FT(0)
        for i in eachindex(new_floes)
            sqr_dist = (centroid[1] - new_floes.centroid[i][1])^2 +
                (centroid[2] -  new_floes.centroid[i][2])^2
            total_inertia += new_floes.moment[i] + new_floes.mass[i] * sqr_dist
            append!(new_floes.parent_id[i], id_lst)
        end

        # Add original floe values, weighted by mass
        for f in eachindex(floes)
            mass_frac = f.mass / total_mass
            moment_frac = f.moment / total_inertia
            new_floes.u .+= f.u * mass_frac
            new_floes.v .+= f.v * mass_frac
            new_floes.ξ .+= f.ξ * moment_frac
            new_floes.p_dudt .+= f.p_dudt * mass_frac
            new_floes.p_dvdt .+= f.p_dvdt * mass_frac
            new_floes.p_dxdt .+= f.p_dxdt * mass_frac
            new_floes.p_dydt .+= f.p_dydt * mass_frac
            new_floes.p_dξdt .+= f.p_dξdt * moment_frac
            # Update stress history
            new_floes.stress_history[1].cb .+= f.stress_history.cb * mass_frac
            new_floes.stress_history[1].total += f.stress_history.total * mass_frac
            for sh in @view new_floes.stress_history[2:end]
                sh.cb .= new_floes.stress_history[1].cb
                sh.total = new_floes.stress_history[1].total
            end
        end
        # new_floes.id ??
        floes.alive .= false  # going to switch this to 'remove' once status is updated
    end
    return new_floes
end