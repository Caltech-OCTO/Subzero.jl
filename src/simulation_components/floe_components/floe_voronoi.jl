"""
    generate_voronoi_coords(
        desired_points,
        scale_fac,
        trans_vec,
        domain_coords,
        rng;
        max_tries = 10,
    )

Generate voronoi coords within a bounding box defined by its lower left corner
and its height and width. Attempt to generate `npieces` cells within the box.
Inputs:
    desired_points  <Int> desired number of voronoi cells
    scale_fac       <Vector{AbstractFloat}> width and height of bounding box -
                        formatted as [w, h] 
    trans_vec       <Vector{AbstractFloat}> lower left corner of bounding box -
                        formatted as [x, y] 
    domain_coords   <Vector{PolyVec{AbstractFloat}}> multipolygon that will
                        eventually be filled with/intersected with the voronoi
                        cells - such as topography
    rng             <RNG> random number generator to generate voronoi cells
    min_to_warn     <Int> minimum number of points to warn if not generated to
                        seed voronoi
    max_tries       <Int> number of tires to generate desired number of points
                        within domain_coords to seed voronoi cell creation
Outputs:
    coords  <Vector{PolyVec{Float}}> vector of polygon coordinates generated by
        voronoi tesselation. These polygons all fall within the space defined by
        the domain_coords. If less polygons than min_to_warn are generated, the
        user will be warned. 
"""
function generate_voronoi_coords(
    desired_points::Int,
    scale_fac,
    trans_vec,
    domain_coords::Vector{<:PolyVec{<:FT}},
    rng,
    min_to_warn::Int;
    max_tries::Int = 10,
) where {FT <: AbstractFloat}
    xpoints = Vector{Float64}()
    ypoints = Vector{Float64}()
    domain_poly = make_multipolygon(GO.tuples(domain_coords))
    area_frac = GO.area(domain_poly) / reduce(*, scale_fac)
    # Increase the number of points based on availible percent of bounding box
    npoints = ceil(Int, desired_points / area_frac)
    current_points = 0
    tries = 0
    while current_points < desired_points && tries <= max_tries
        x = rand(rng, npoints)
        y = rand(rng, npoints)
        # Check which of the scaled and translated points are within the domain coords
        in_idx = [GO.coveredby(
            (scale_fac[1] * x[i] .+ trans_vec[1], scale_fac[2] * y[i] .+ trans_vec[2]),
            domain_poly
        ) for i in eachindex(x)]
        current_points += sum(in_idx)
        tries += 1
        append!(xpoints, x[in_idx])
        append!(ypoints, y[in_idx])
    end
    # If we generated too many cells, remove extra
    if current_points > desired_points
        xpoints = xpoints[1:desired_points]
        ypoints = ypoints[1:desired_points]
        current_points = desired_points
    end
    # Warn if we didn't generate enough cells
    if current_points < min_to_warn
        @warn "Only $current_points floes were able to be generated in \
            $max_tries tries during voronoi tesselation."
    end
    # Make voronoi cells into floes
    if current_points > 1
        tess_cells = voronoicells(
            xpoints,
            ypoints,
            Rectangle(GB.Point2((0.0, 0.0)), GB.Point2((1.0, 1.0))),
            rng = rng
        ).Cells
        # Scale and translate voronoi coordinates
        tcoords = Vector{PolyVec{FT}}(undef, length(tess_cells))
        for i in eachindex(tess_cells)
            tcoords[i] = [valid_ringvec!([
                Vector(c) .* scale_fac .+ trans_vec
                for c in tess_cells[i]
            ])]
        end
        return tcoords
    else
        return Vector{PolyVec{FT}}()
    end
end

function _rand_points_in_poly(rng, poly, npoints, max_tries)
    current_points = 0
    tries = 0

end