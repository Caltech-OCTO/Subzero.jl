using GeoInterface
using MAT
using JLD2

function polygons_to_matlab(jl_file::Vector, mat_file::String)
    # get julia polygons where jl_file is a floes file output by a simulation
    f = jldopen(jl_file)
    polys = f["poly"]["0"]
    close(f)

    n = length(polys)
    # MATLAB-style cell array: Vector of matrices
    poly_cells = Vector{Matrix{Float64}}(undef, n)

    for i in 1:n
        poly = polys[i]

        # Get exterior ring
        ring = GeoInterface.getexterior(poly)

        # Extract coordinates
        coords = collect(GeoInterface.getpoint.(ring))

        # Split into x/y
        x = [c[1] for c in coords]
        y = [c[2] for c in coords]

        # MATLAB expects Nx2 for polyshape
        poly_cells[i] = hcat(x, y)
    end

    # Write MAT file
    matwrite(mat_file, Dict("polygons" => poly_cells))
end
