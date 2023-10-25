
function bin_floe_centroids(floes, grid, Nx, Ny)
    nfloes = length(floes)
    Matrix{Vector{Int}}(zeros(Int, nfloes / Nx*Ny), Nx, Ny)
    floe_bins = zeros(Int, Nx, Ny, nfloes / Nx*Ny)
    Lx = grid.xf - grid.x0
    Ly = grid.yf - grid.y0
    Δx = Lx / Nx
    Δy = Ly / Ny
    for i in eachindex(floes)
        xp, yp = floes.centroid[i]
        xidx = floor(Int, (xp - grid.x0) / Δx) + 1
        yidx = floor(Int, (yp - grid.y0) / Δy) + 1

    end
end


function timestep_welding!(floes)

    return 
end