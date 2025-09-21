export make_polygon

"""
    valid_ringvec(coords::RingVec{FT})

Takes a RingVec object and make sure that the last element has the same first
element as last element and that other than these two elements there are no
duplicate, adjacent vertices. Also asserts that the ring as at least three
elements or else it cannot be made into a valid ring as it is a line segment. 
"""
function valid_ringvec!(ring)
    deleteat!(ring, findall(i->ring[i]==ring[i+1], collect(1:length(ring)-1)))
    if ring[1] != ring[end]
        push!(ring, deepcopy(ring[1]))
    end
    @assert length(ring) > 3 "Polgon needs at least 3 distinct points."
    return ring
end

"""
    valid_polyvec(coords::PolyVec{FT})

Takes a PolyVec object and make sure that the last element of each "ring"
(vector of vector of floats) has the same first element as last element and has
not duplicate adjacent elements. Also asserts that each "ring" as at least three
distinct elements or else it is not a valid ring, but rather a line segment. 
"""
function valid_polyvec!(coords)
    for ring in coords
        valid_ringvec!(ring)
    end
    return coords
end

"""
    intersect_polys(p1, p2)

Intersect two geometries and return a list of polygons resulting.
Inputs:
    p1  <AbstractGeometry>
    p2  <AbstractGeometry>
Output:
    Vector of Polygons
"""
intersect_polys(p1, p2, ::Type{FT} = Float64; kwargs...) where FT = GO.intersection(p1, p2, FT; target = GI.PolygonTrait(), fix_multipoly = nothing)
diff_polys(p1, p2, ::Type{FT} = Float64; kwargs...) where FT = GO.difference(p1, p2, FT; target = GI.PolygonTrait(), fix_multipoly = nothing) 
union_polys(p1, p2, ::Type{FT} = Float64; kwargs...) where FT = GO.union(p1, p2, FT; target = GI.PolygonTrait(), fix_multipoly = nothing)
simplify_poly(p, tol) = GO.simplify(p; tol = tol)

function _translate_poly(::Type{FT}, p, Δx, Δy) where FT
    t = CoordinateTransformations.Translation(Δx, Δy)
    # TODO: can remove the tuples call after GO SVPoint PR
    return GO.tuples(GO.transform(t, p), FT)
end

function _move_poly(::Type{FT}, poly, Δx, Δy, Δα, cx = zero(FT), cy = zero(FT)) where FT
    rot = CoordinateTransformations.LinearMap(Rotations.Angle2d(Δα))
    cent_rot = CoordinateTransformations.recenter(rot, (cx, cy))
    trans = CoordinateTransformations.Translation(Δx, Δy)
    # TODO: can remove the tuples call after GO SVPoint PR
    return GO.tuples(GO.transform(trans ∘ cent_rot, poly), FT)::Polys{FT}
end

make_polygon(coords::PolyVec) = GI.Polygon(GO.tuples(coords))
make_polygon(tuple_coords) = GI.Polygon(tuple_coords)
make_polygon(ring::GI.LinearRing) = GI.Polygon([ring])
make_multipolygon(coords::Vector{<:PolyVec}) = GI.MultiPolygon(GO.tuples(coords))
make_multipolygon(tuple_coords) = GI.MultiPolygon(tuple_coords)
make_multipolygon(polys::Vector{<:GI.Polygon}) = GI.MultiPolygon(polys)
function make_multipolygon(polys::Vector{<:StaticQuadrilateral{FT}}) where FT
    new_polys = Vector{Polys{FT}}(undef, length(polys))
    for (i, poly) in enumerate(polys)
        new_polys[i] = make_polygon([[p for p in GI.getpoint(poly)]])
    end
    return make_multipolygon(new_polys)
end

function _make_bounding_box_polygon(::Type{FT}, xmin, xmax, ymin, ymax) where FT
    points = ((xmin, ymin),  (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin))
    ring = GI.LinearRing(SA.SVector{5, Tuple{FT, FT}}(points))
    return  GI.Polygon(SA.SVector(ring))
end

"""
    hashole(poly::Polys)

Determine if polygon has one or more holes
Inputs:
    poly <Polygon> polygon
Outputs:
    <Bool> true if there is a hole in the polygons, else false
"""
function hashole(poly::Polys)
    return GI.nhole(poly) > 0
end 

function rmholes!(poly::Polys)
    deleteat!(poly.geom, 2:GI.nring(poly))
    return
end

# Find the length of the maximum radius of a given polygon
function _calc_max_radius(poly, cent, ::Type{T}) where T
    max_rad_sqrd = zero(T)
    Δx, Δy = GO._tuple_point(cent, T)
    for pt in GI.getpoint(GI.getexterior(poly))
        x, y = GO._tuple_point(pt, T)
        x, y = x - Δx, y - Δy
        rad_sqrd = x^2 + y^2
        if rad_sqrd > max_rad_sqrd
            max_rad_sqrd = rad_sqrd
        end
    end
    return sqrt(max_rad_sqrd)
end

"""
    which_vertices_match_points(ipoints, coords, atol)

Find which vertices in coords match given points
Inputs:
    points <Vector{Tuple{Float, Float} or Vector{Vector{Float}}}> points to
                match to vertices within polygon
    region  <Polygon> polygon 
    atol    <Float> distance vertex can be away from target point before being
                classified as different points
Output:
    Vector{Int} indices of points in polygon that match the intersection points
Note: 
    If last coordinate is a repeat of first coordinate, last coordinate index is
    NOT recorded.
"""
function which_vertices_match_points(points, region::Polys{FT}, atol = 1) where FT
    idxs = Vector{Int}()
    npoints = length(points)
    if points[1] == points[end]
        npoints -= 1
    end
    for i in 1:npoints  # find which vertex matches point
        min_dist = FT(Inf)
        min_vert = 1
        for (j, pt) in enumerate(GI.getpoint(GI.getexterior(region)))
            dist = sqrt(GO.distance(pt, points[i], FT))
            if dist < min_dist
                min_dist = dist
                min_vert = j
            end
        end
        if min_dist < atol
            push!(idxs, min_vert)
        end
    end
    return sort!(idxs)
end

#=
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
=#
function _generate_voronoi_coords(  # TODO: maybe move to floe utils since it is used in mutliple places!
    ::Type{FT},
    desired_points::Int,
    scale_fac,
    trans_vec,
    domain_poly,
    rng,
    min_to_warn::Int;
    max_tries::Int = 10,
) where {FT <: AbstractFloat}
    xpoints = Vector{Float64}()
    ypoints = Vector{Float64}()
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