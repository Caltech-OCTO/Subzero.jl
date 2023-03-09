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