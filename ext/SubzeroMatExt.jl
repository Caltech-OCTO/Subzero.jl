module SubzeroMatExt

using Subzero, MAT, JLD2, GeoInterface
import Subzero: polygons_to_matlab

function polygons_to_matlab(jl_file, mat_file)
    # get julia polygons where jl_file is a floes file output by a simulation
    f = jldopen(jl_file)
    polys = f["poly"]["0"]
    close(f)

    n = length(polys)
    # MATLAB-style cell array: Vector of matrices
    poly_cells = Vector{Matrix{Float64}}(undef, n)

    for i in range(1, n)
        poly = polys[i]

        # Get exterior ring
        ring = GeoInterface.getexterior(poly)

        # Extract coordinates
        coords = GeoInterface.getpoint(ring)

        # Split into x/y
        x = first.(coords)
        y = last.(coords)

        # MATLAB expects Nx2 for polyshape
        poly_cells[i] = hcat(x, y)
    end

    # Write MAT file
    matwrite(mat_file, Dict("polygons" => poly_cells))
end

end  # module