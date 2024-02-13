
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
    find_poly_centroid(poly)

Syntactic sugar for using LibGEOS to find a polygon's centroid
Input:
    poly    <LibGEOS.Polygon or LibGEOS. MultiPolygon>
Output:
    Vector{Float64} [x, y] where this represents the centroid in the xy plane
"""
find_poly_centroid(poly) =
    LG.GeoInterface.coordinates(LG.centroid(poly))::Vector{Float64}

"""
    find_poly_coords(poly)

Syntactic sugar for using LibGEOS to find a polygon's coordinates
Input:
    poly    <LibGEOS.Polygon>
Output:
    <PolyVec> representing the floe's coordinates xy plane
"""
find_poly_coords(poly::Union{LG.Polygon, GI.Polygon}) =
    LG.GeoInterface.coordinates(poly)

"""
    get_polygons(geom)

Returns an empty polygon list as non-polygon element was provided
Inputs:
    geom    <LG.AbstractGeometry>
Outputs:
    <Vector{LG.Polygon}>
"""
function get_polygons(geom, ::Type{T} = Float64) where T
    polys =  Vector{Polys{T}}()
    _get_polygons!(geom, polys)
    return polys
end

_get_polygons!(_, polys) = nothing

function _get_polygons!(
    geom::Polys{T},
    polys,
) where T
    if !GI.isempty(geom) && GO.area(geom) > 0
        push!(polys, geom)
    end
    return
end

function _get_polygons!(
    collection::Union{MultiPolys{T}, LG.GeometryCollection},
    polys,
) where T
    for geom in GI.getgeom(collection)
        _get_polygons!(geom, polys)
    end
    return
end

"""
    intersect_polys(p1, p2)

Intersect two geometries and return a list of polygons resulting.
Inputs:
    p1  <LG.AbstractGeometry>
    p2  <LG.AbstractGeometry>
Output:
    Vector of LibGEOS Polygons
"""
intersect_polys(p1, p2) = get_polygons(
    LG.intersection(
        p1,
        p2,
    )::LG.Geometry,
)

"""
    intersect_coords(c1, c2)

Intersect geometries using their coordinates to return list of resulting
polygons.
Inputs:
    c1  <PolyVec>
    c2  <PolyVec>
Output:
    Vector of LibGEOS Polygons
"""
intersect_coords(c1, c2) = intersect_polys(
    LG.Polygon(c1),
    LG.Polygon(c2),
)


"""
    polyvec_extrema(coords)

Finds extremal x and y values for given coordiantes, defining a tight bounding
box. 
Inputs:
    coords <PolyVec>
Outputs: xmin, xmax, ymin, ymax
"""
function polyvec_extrema(coords::PolyVec{FT}) where FT
    xmin = ymin = FT(Inf)
    xmax = ymax = FT(-Inf)
    for i in eachindex(coords[1])
        x, y = coords[1][i]
        if x < xmin
            xmin = x
        end
        if x > xmax
            xmax = x
        end
        if y < ymin
            ymin = y
        end
        if y > ymax
            ymax = y
        end
    end
    return xmin, xmax, ymin, ymax
end

"""
    deepcopy_floe(floe::LazyRow{Floe{FT}})

Deepcopy of a floe by creating a new floe and copying all fields.
Inputs:
    floe    <Floe>
Outputs:
    New floe with floes that are equal in value. Any vector fields are copies so
    they share values, but not referance.
"""
function deepcopy_floe(floe::LazyRow{Floe{FT}}) where {FT}
    f = Floe{FT}(
        centroid = copy(floe.centroid),
        coords = translate(floe.coords, 0, 0),
        height = floe.height,
        area = floe.area,
        mass = floe.mass,
        rmax = floe.rmax,
        moment = floe.moment,
        angles = copy(floe.angles),
        x_subfloe_points = copy(floe.x_subfloe_points),
        y_subfloe_points = copy(floe.y_subfloe_points),
        α = floe.α,
        u = floe.u,
        v = floe.v,
        ξ = floe.ξ,
        status = Status(floe.status.tag, copy(floe.status.fuse_idx)),
        id = floe.id,
        ghost_id = floe.ghost_id,
        parent_ids = copy(floe.parent_ids),
        ghosts = copy(floe.ghosts),
        fxOA= floe.fxOA,
        fyOA = floe.fyOA,
        trqOA = floe.trqOA,
        hflx_factor = floe.hflx_factor,
        overarea = floe.overarea,
        collision_force = copy(floe.collision_force),
        collision_trq = floe.collision_trq,
        stress = copy(floe.stress),
        stress_history = StressCircularBuffer{FT}(
            capacity(floe.stress_history.cb),
        ),
        strain = copy(floe.strain),
        p_dxdt = floe.p_dxdt,
        p_dydt = floe.p_dydt,
        p_dudt = floe.p_dudt,
        p_dvdt = floe.p_dvdt,
        p_dξdt = floe.p_dξdt,
        p_dαdt = floe.p_dαdt,
    )
    f.stress_history.total .= copy(floe.stress_history.total)
    append!(f.stress_history.cb, copy(floe.stress_history.cb))
    return f
end

"""
    find_multipoly_coords(poly)

Syntactic sugar for using LibGEOS to find a multipolygon's coordinates
Input:
    poly    <LibGEOS.Polygon>
Output:
    <Vector{PolyVec}> representing the floe's coordinates xy plane
"""
find_multipoly_coords(poly::Union{LG.Polygon, GI.Polygon}) =
    [find_poly_coords(poly)]


"""
    find_multipoly_coords(multipoly)

Syntactic sugar for using LibGEOS to find a multi-polygon's coordinates
Input:
    poly    <LibGEOS. MultiPolygon>
Output:
    <Vector{PolyVec}> representing the floe's coordinates xy plane
"""
find_multipoly_coords(multipoly::LG.MultiPolygon) =
    LG.GeoInterface.coordinates(multipoly)::Vector{PolyVec{Float64}}

"""
    translate!(coords, Δx, Δy)

Make a copy of given coordinates and translate by given deltas. 
Inputs:
    coords PolyVec{Float}
    vec <Vector{Real}>
Output:
    Updates given coords
"""
function translate(coords::PolyVec{FT}, Δx, Δy) where {FT<:AbstractFloat}
    new_coords = [[Vector{Float64}(undef, 2) for _ in eachindex(coords[1])]]
    for i in eachindex(coords[1])
        new_coords[1][i][1] = coords[1][i][1] + Δx
        new_coords[1][i][2] = coords[1][i][2] + Δy 
    end
    return new_coords
end

"""
    translate!(coords, Δx, Δy)

Translate each of the given coodinates by given deltas in place
Inputs:
    coords PolyVec{Float}
    vec <Vector{Real}>
Output:
    Updates given coords
"""
function translate!(coords::PolyVec{FT}, Δx, Δy) where {FT<:AbstractFloat}
    for i in eachindex(coords)
        for j in eachindex(coords[i])
            coords[i][j][1] += Δx
            coords[i][j][2] += Δy
        end
    end
    return
end

"""
    rotate_radians!(coords::PolyVec, α)

Rotate a polygon's coordinates by α radians around the origin.
Inputs:
    coords  <PolyVec{AbstractFloat}> polygon coordinates
    α       <Real> radians to rotate the coordinates
Outputs:
    Updates coordinates in place
"""
function rotate_radians!(coords::PolyVec, α)
    for i in eachindex(coords)
        for j in eachindex(coords[i])
            x, y = coords[i][j]
            coords[i][j][1] = cos(α)*x - sin(α)*y
            coords[i][j][2] = sin(α)*x + cos(α)*y
        end
    end
    return
end

"""
    rotate_degrees!(coords::PolyVec, α)

Rotate a polygon's coordinates by α degrees around the origin.
Inputs:
    coords <PolyVec{AbstractFloat}> polygon coordinates
    α       <Real> degrees to rotate the coordinates
Outputs:
    Updates coordinates in place
"""
rotate_degrees!(coords::PolyVec, α) = rotate_radians!(coords, α * π/180)

"""
    hashole(coords::PolyVec{FT})

Determine if polygon coordinates have one or more holes
Inputs:
    coords <PolyVec{Float}>
Outputs:
    <Bool>
"""
function hashole(coords::PolyVec{FT}) where FT<:AbstractFloat
    return length(coords) > 1
end

"""
    hashole(poly::LG.Polygon)

Determine if polygon has one or more holes
Inputs:
    poly <LibGEOS.Polygon> LibGEOS polygon
Outputs:
    <Bool> true if there is a hole in the polygons, else false
"""
function hashole(poly::Union{LG.Polygon, GI.Polygon})
    return GI.nhole(poly) > 0
end 

"""
    hashole(multipoly::LG.MultiPolygon)

Determine if any of multipolygon's internal polygons has holes
Inputs:
    multipoly <LibGEOS.MultiPolygon> LibGEOS multipolygon
Outputs:
    <Bool> true if there is a hole in any of the polygons, else false
"""
function hashole(multipoly::LG.MultiPolygon)
    poly_lst = LG.getGeometries(multipoly)::Vector{LG.Polygon}
    for poly in poly_lst
        if hashole(poly)
            return true
        end
    end
    return false
end

"""
    rmholes(coords::PolyVec{FT})

Remove polygon coordinates's holes if they exist
Inputs:
    coords <PolyVec{Float}> polygon coordinates
Outputs:
    <PolyVec{Float}> polygonc coordinates after removing holes
"""
function rmholes(coords::PolyVec{FT}) where {FT<:AbstractFloat}
    return [coords[1]]
end

function rmholes!(coords::PolyVec{FT}) where {FT<:AbstractFloat}
    if length(coords) > 1
        deleteat!(coords, 2:length(coords))
    end
end

"""
    rmholes(poly::LG.Polygon)

Remove polygon's holes if they exist
Inputs:
    poly <LibGEOS.Polygon> LibGEOS polygon
Outputs:
    <LibGEOS.Polygon>  LibGEOS polygon without any holes
"""
function rmholes(poly::Union{LG.Polygon, GI.Polygon})
    if hashole(poly)
        return LG.Polygon(LG.exteriorRing(poly))
    end
    return poly
end

"""
    rmholes(multipoly::LG.MultiPolygon)

Remove holes from each polygon of a multipolygon if they exist
Inputs:
    multipoly <LibGEOS.MultiPolygon> multipolygon coordinates
Outputs:
    <LibGEOS.MultiPolygon> multipolygon without any holes
"""
function rmholes(multipoly::LG.MultiPolygon)
    poly_lst = LG.getGeometries(multipoly)::Vector{LG.Polygon}
    nohole_lst = LG.Polygon[]
    for poly in poly_lst
        push!(nohole_lst, rmholes(poly))
    end
    return LG.MultiPolygon(nohole_lst)
end

"""
    sortregions(poly::LG.Polygon)

Returns given polygon within a vector as it is the only region
Inputs:
    poly <LibGEOS.Polygon> LibGEOS polygon
Outputs:
    <Vector{LibGEOS.Polygon}> single LibGEOS polygon within list
"""
function sortregions(poly::Union{LG.Polygon, GI.Polygon})
    return [poly]
end

"""
    sortregions(multipoly::LG.MultiPolygon)

Sorts polygons within a multi-polygon by area in descending order
Inputs:
    multipoly <LibGEOS.MultiPolygon> LibGEOS multipolygon
Outputs:
    <Vector{LibGEOS.Polygon}> list of LibGEOS polygons sorted in descending
        order by area
"""
function sortregions(multipoly::LG.MultiPolygon)
    return sort!(LG.getGeometries(multipoly), by=GO.area, rev=true)
end

"""
    separate_xy(coords::PolyVec{T})

Pulls x and y coordinates from standard polygon vector coordinates into seperate
vectors. Only keeps external coordinates and disregards holes
Inputs:
    coords <PolyVec{Float}> polygon coordinates
Outputs: 
    x <Vector{Float}> x coordinates
    y <Vector{Float}> y coordinates
"""
function separate_xy(coords::PolyVec{T}) where {T<:AbstractFloat}
    x = first.(coords[1])
    y = last.(coords[1])
    return x, y 
end

"""
    calc_moment_inertia(coords::PolyVec{T}, h; rhoice = 920.0)

Calculate the mass moment of intertia from polygon coordinates.
Inputs:
    coords      <PolyVec{Float}>
    h           <Real> height of floe
    rhoice      <Real> Density of ice
Output:
    <Float> mass moment of inertia
Note:
    Assumes that first and last point within the coordinates are the same.
    Will not give correct answer otherwise.

Based on paper: Marin, Joaquin."Computing columns, footings and gates through
moments of area." Computers & Structures 18.2 (1984): 343-349.
"""
function calc_moment_inertia(
    coords::PolyVec{T},
    centroid,
    h;
    ρi = 920.0
) where {T<:AbstractFloat}
    x, y = separate_xy(coords)
    x .-= centroid[1]
    y .-= centroid[2]
    N = length(x)
    wi = x[1:N-1] .* y[2:N] - x[2:N] .* y[1:N-1]
    Ixx = 1/12 * sum(wi .* ((y[1:N-1] + y[2:N]).^2 - y[1:N-1] .* y[2:N]))
    Iyy = 1/12 * sum(wi .* ((x[1:N-1] + x[2:N]).^2 - x[1:N-1] .* x[2:N]))
    return abs(Ixx + Iyy)*h*ρi;
end

"""
    calc_moment_inertia(poly::LG.Polygon, h; rhoice = 920.0)

Calculate the mass moment of intertia from a LibGEOS polygon object using above
coordinate-based moment of intertia function.
Inputs:
    poly      LibGEOS.Polygon
    h           <Real> height of floe
    rhoice      <Real> Density of ice
Output:
    <Float> mass moment of inertia
"""
calc_moment_inertia(poly::Union{LG.Polygon, GI.Polygon}, h; ρi = 920.0) = 
    calc_moment_inertia(
        find_poly_coords(poly),
        find_poly_centroid(poly),
        h,
        ρi = ρi,
    )

"""
    orient_coords(coords)

Take given coordinates and make it so that the first point has the smallest
x-coordiante and so that the coordinates are ordered in a clockwise sequence.
Duplicates vertices will be removed and the coordiantes will be closed (first
and last point are the same).

Input:
    coords  <RingVec> vector of points [x, y]
Output:
    coords  <RingVec> oriented clockwise with smallest x-coordinate first
"""
function orient_coords(coords::RingVec)
    # extreem_idx is point with smallest x-value - if tie, choose lowest y-value
    extreem_idx = 1
    for i in eachindex(coords)
        ipoint = coords[i]
        epoint = coords[extreem_idx]
        if ipoint[1] < epoint[1]
            extreem_idx = i
        elseif ipoint[1] == epoint[1] && ipoint[2] < epoint[2]
            extreem_idx = i
        end
    end
    # extreem point must be first point in list
    new_coords = similar(coords)
    circshift!(new_coords, coords, -extreem_idx + 1)
    valid_ringvec!(new_coords)

    # if coords are counterclockwise, switch to clockwise
    orient_matrix = hcat(
        ones(3),
        vcat(new_coords[1]', new_coords[2]', new_coords[end-1]') # extreem/adjacent points
    )
    if det(orient_matrix) > 0
        reverse!(new_coords)
    end
    return new_coords
end

"""
    calc_point_poly_dist(xp::Vector{T},yp::Vector{T}, vec_poly::PolyVec{T})

Compute the distances from each one of a set of np points on a 2D plane to a
polygon. Distance from point j to an edge k is defined as a distance from this
point to a straight line passing through vertices v(k) and v(k+1), when the
projection of point j on this line falls INSIDE segment k; and to the closest of
v(k) or v(k+1) vertices, when the projection falls OUTSIDE segment k.
Inputs:
    xp  <Vector{Float}> x-coordinates of points to find distance from vec_poly
    yp  <Vector{Float}> y-coordiantes of points to find distance from vec_poly
    vec_poly    <PolyVec{Float}> coordinates of polygon
Outputs:
    <Vector{AbstractFloat}>List of distances from each point to the polygon. If
    the point is inside of the polygon the value will be negative. This does not
    take holes into consideration.

Note - Translated into Julia from the following program:
p_poly_dist by Michael Yoshpe - last updated in 2006.
We mimic version 1 functionality with 4 inputs and 1 output.
Only needed code was translated.
"""
function calc_point_poly_dist(
    xp::Vector{FT},
    yp::Vector{FT},
    vec_poly::PolyVec{FT}
) where {FT<:AbstractFloat}
    @assert length(xp) == length(yp)
    min_lst = if !isempty(xp)
        # Vertices in polygon and given points
        Pv = reduce(hcat, valid_polyvec!(vec_poly)[1])'
        Pp = hcat(xp, yp)
        np = length(xp)
        nv = length(vec_poly[1])
        # Distances between all points and vertices in x and y
        x_dist = repeat(Pv[:, 1], 1, np)' .- repeat(Pp[:, 1], 1, nv)
        y_dist = repeat(Pv[:, 2], 1, np)' .- repeat(Pp[:, 2], 1, nv)
        p2c_dist = hypot.(x_dist, y_dist)
        # minimum distance to vertices
        min_dist, min_idx = findmin(p2c_dist, dims = 2)
        # Coordinates of consecutive vertices
        V1 = Pv[1:end-1, :]
        V2 = Pv[2:end, :]
        Δv = V2 .-  V1
        # Vector of distances between each pair of consecutive vertices
        vds = hypot.(Δv[:, 1], Δv[:, 2])

        if (cumsum(vds)[end-1] - vds[end]) < 10eps(FT)
            throw(ArgumentError("Polygon vertices should not lie on a straight \
                line"))
        end

        #= Each pair of consecutive vertices V1[j], V2[j] defines a rotated
        coordinate system with origin at V1[j], and x axis along the vector
        V2[j]-V1[j]. cθ and sθ rotate from original to rotated system =#
        cθ = Δv[:, 1] ./ vds
        sθ = Δv[:, 2] ./  vds
        Cer = zeros(FT, 2, 2, nv-1)
        Cer[1, 1, :] .= cθ
        Cer[1, 2, :] .= sθ
        Cer[2, 1, :] .= -sθ
        Cer[2, 2, :] .= cθ

        # Build origin translation vector P1r in rotated frame by rotating V1
        V1r = hcat(cθ .* V1[:, 1] .+ sθ .* V1[:, 2], 
            -sθ .* V1[:, 1] .+ cθ .* V1[:, 2])

        #= Ppr is a 3D array of size 2*np*(nv-1). Ppr(1,j,k) is an X coordinate
        of point j in coordinate systems defined by segment k. Ppr(2,j,k) is its
        Y coordinate. =#
        Ppr = zeros(FT, 2, np, nv-1)
        # Rotation and Translation
        Ppr[1, :, :] .= Pp * Cer[1, :, :] .-
            permutedims(repeat(V1r[:, 1], 1, 1, np), [2, 3, 1])[1, :, :]
        Ppr[2, :, :] .= Pp * Cer[2, :, :] .-
            permutedims(repeat(V1r[:, 2], 1, 1, np), [2, 3, 1])[1, :, :]

        # x and y coordinates of the projected (cross-over) points in original
        # coordinate
        r = Ppr[1, :, :]
        cr = Ppr[2, :, :]
        B = fill(convert(FT, Inf), np, nv-1)
        #= For the projections that fall inside the segments, find the minimum
        distances from points to their projections (note, that for some points
        these might not exist) =#
        for i in eachindex(r)
            if r[i] > 0 && r[i] < vds[cld(i, np)]
                B[i] = cr[i]
            end
        end
        cr_min, cr_min_idx = findmin(abs.(B), dims = 2)
        #= For projections that fall outside segments, closest point is a vertex
        These points have a negative value if point is actually outside of
        polygon =#
        in_poly = inpoly2(Pp, Pv)
        dmin = cr_min
        for i in eachindex(dmin)
            if isinf(dmin[i]) || (cr_min_idx[i] != min_idx[i] && cr_min[i] > min_dist[i])
                dmin[i] = min_dist[i]
            end
            if in_poly[i, 1] ||  in_poly[i, 2]
                dmin[i] *= -1
            end
        end
        dmin[:, 1]  # Turn array into a vector
    else
        FT[]
    end
    return min_lst
end

"""
    intersect_lines(l1, l2)

Finds the intersection points of two curves l1 and l2. The curves l1, l2 can be
either closed or open. In this version, l1 and l2 must be distinct. If no
intersections are found, the returned P is empty.
Inputs:
    l1 <PolyVec{Float}> line/polygon coordinates
    l2 <PolyVec{Float}> line/polygon coordinates
Outputs:
    <Set{Tuple{Float, Float}}> Set of points that are at the intersection of the
        two line segments.
"""
function intersect_lines(l1::PolyVec{FT}, l2) where {FT}
    points = Vector{Tuple{FT, FT}}()
    for i in 1:(length(l1[1]) - 1)
         # First line represented as a1x + b1y = c1
         x11 = l1[1][i][1]
         x12 = l1[1][i+1][1]
         y11 = l1[1][i][2]
         y12 = l1[1][i+1][2]
         a1 = y12 - y11
         b1 = x11 - x12
         c1 = a1 * x11 + b1 * y11
        for j in 1:(length(l2[1]) - 1)
            # Second line represented as a2x + b2y = c2
            x21 = l2[1][j][1]
            x22 = l2[1][j+1][1]
            y21 = l2[1][j][2]
            y22 = l2[1][j+1][2]
            a2 = y22 - y21
            b2 = x21 - x22
            c2 = a2 * x21 + b2 * y21

            determinant = a1 * b2 - a2 * b1
            # Find place there two lines cross
            # Note that lines extend beyond given line segments
            if determinant != 0
                x = (b2*c1 - b1*c2)/determinant
                y = (a1*c2 - a2*c1)/determinant
                p = (x, y)
                # Make sure intersection is on given line segments
                if min(x11, x12) <= x <= max(x11, x12) &&
                    min(x21, x22) <= x <= max(x21, x22) &&
                    min(y11, y12) <= y <= max(y11, y12) &&
                    min(y21, y22) <= y <= max(y21, y22) &&
                    !(p in points)
                    push!(points, p)
                end
            end
        end
    end
    return points
end

"""
    which_vertices_match_points(ipoints, coords, atol)

Find which vertices in coords match given points
Inputs:
    points <Vector{Tuple{Float, Float} or Vector{Vector{Float}}}> points to
                match to vertices within polygon
    coords  <PolVec> polygon coordinates
    atol    <Float> distance vertex can be away from target point before being
                classified as different points
Output:
    Vector{Int} indices of points in polygon that match the intersection points
Note: 
    If last coordinate is a repeat of first coordinate, last coordinate index is
    NOT recorded.
"""
function which_vertices_match_points(
    points,
    coords::PolyVec{FT},
    atol = 1,
) where {FT}
    idxs = Vector{Int}()
    npoints = length(points)
    if points[1] == points[end]
        npoints -= 1
    end
    @views for i in 1:npoints  # find which vertex matches point
        min_dist = FT(Inf)
        min_vert = 1
        for j in eachindex(coords[1])
            dist = sqrt(
                (coords[1][j][1] - points[i][1])^2 +
                (coords[1][j][2] - points[i][2])^2,
            )
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

"""
euclidian_dist(c, idx2, idx1)

Calculate euclidean distance between two points within given coordinates
"""
euclidian_dist(c, idx2, idx1) = sqrt(
    (c[1][idx2][1] - c[1][idx1][1])^2 +
    (c[1][idx2][2] - c[1][idx1][2])^2 
)

"""
    which_points_on_edges(points, coords; atol = 1e-1)

Find which points are on the coordinates of the given polygon.
Inputs:
    points <Vector{Tuple{Float, Float} or Vector{Vector{Float}}}> points to
        match to edges within polygon
    coords  <PolVec> polygon coordinates
    atol    <Float> distance target point can be from an edge before being
                classified as not on the edge
"""
function which_points_on_edges(points, coords; atol = 1e-1)
    idxs = Vector{Int}()
    nedges = length(coords[1]) - 1
    npoints = length(points)
    if points[1] == points[end]
        npoints -= 1
    end
    for i in 1:nedges
        x1, y1 = coords[1][i]
        x2, y2 = coords[1][i+1]
        Δx_edge = x2 - x1
        Δy_edge = y2 - y1
        for j in 1:npoints
            xp, yp = points[j]
            Δx_point = xp - x1
            Δy_point = yp - y1
            if (!(j in idxs) && (
                (  # vertical edge
                    isapprox(Δx_edge, 0, atol = atol) &&
                    isapprox(Δx_point, 0, atol = atol) &&
                    0 < Δy_point / Δy_edge < 1
                ) ||
                (  # horizontal edge
                    isapprox(Δy_edge, 0, atol = atol) &&
                    isapprox(Δy_point, 0, atol = atol)
                    && 0 < Δx_point / Δx_edge < 1) ||
                (  # point is a vertex
                    isapprox(Δx_point, 0, atol = atol) &&
                    isapprox(Δy_point, 0, atol = atol)
                ) ||
                ( # has the same slope and is between edge points
                    isapprox(Δy_edge/Δx_edge, Δy_point/Δx_point, atol = atol) &&
                    0 < Δx_point / Δx_edge < 1 && 0 < Δy_point / Δy_edge < 1
                )
            ))
                push!(idxs, j)
            end
        end
    end
    sort!(idxs)
    return idxs
end

"""
    check_for_edge_mid(c, start, stop, shared_idx, shared_dist, running_dist)

Check if indices from start to stop index of given coords includes midpoint
given the shared distance and return midpoint if it exists in given range.
Inputs:
    c               <PolyVec> floe coordinates
    start           <Int> index of shared_index list to start search from
    stop            <Int> index of shared_index list to stop search at
    shared_idx      <Vector{Int}> list of indices of c used to calculate midpoint
    shared_dist     <Float> total length of edges considered from shared_idx
    running_dist    <Float> total length of edges traveled along in midpoint
                        search so far
Outputs:
    mid_x           <Float> x-coordinate of midpoint, Inf if midpoint not in
                        given range
    mid_y           <Float> y-coordinate of midpoint, Inf if midpoint not in
                        given range
    running_dist    <Float> sum of distances travelled along shared edges
"""
function check_for_edge_mid(c, start, stop, shared_idx, shared_dist,
    running_dist::FT,
) where FT
    mid_x = FT(Inf)
    mid_y = FT(Inf)
    for i in start:(stop-1)
        # Indices of c that are endpoints of current edge
        idx1 = shared_idx[i]
        idx2 = shared_idx[i + 1]
        # Lenght of edge
        edge_dist = euclidian_dist(c, idx2, idx1)
        # if midpoint is on current edge
        if running_dist + edge_dist >= shared_dist / 2
            frac = ((shared_dist / 2) - running_dist) / edge_dist
            mid_x = c[1][idx1][1] + (c[1][idx2][1] - c[1][idx1][1]) * frac
            mid_y = c[1][idx1][2] + (c[1][idx2][2] - c[1][idx1][2]) * frac
            break
        else  # move on to the next edge
            running_dist += edge_dist
        end
    end
    return mid_x, mid_y, running_dist
end

"""
    find_shared_edges_midpoint(c1, c2)

Find "midpoint" of shared polygon edges by distance
Inputs:
    c1      <PolVec> polygon coordinates for floe 1
    c2      <PolVec> polygon coordinates for floe 2
Outputs:
    mid_x   <Float> x-coordinate of midpoint
    mid_y   <Float> y-coordinate of midpoint
"""
function find_shared_edges_midpoint(c1::PolyVec{FT}, c2; atol = 1e-1) where {FT}
    # Find which points of c1 are on edges of c2
    shared_idx = which_points_on_edges(
        c1[1],
        c2;
        atol,
    )
    if shared_idx[1] == 1
         # due to repeated first point/last point
        push!(shared_idx, length(c1[1]))
    end
    shared_dist = FT(0)
    gap_idx = 1
    # Determine total length of edges that are shared between the two floes
    nshared_points = length(shared_idx)
    for i in 1:(nshared_points - 1)
        idx1 = shared_idx[i]
        idx2 = shared_idx[i + 1]
        if idx2 - idx1 == 1
            shared_dist += euclidian_dist(c1, idx2, idx1)
        elseif shared_dist > 0
            gap_idx = i + 1
            # Add distance wrapping around from last index to first
            shared_dist += euclidian_dist(c1, shared_idx[end], shared_idx[1])
        end
    end
    # Determine mid-point of shared edges by distance
    running_dist = FT(0)
    mid_x, mid_y, running_dist = check_for_edge_mid(c1, gap_idx, nshared_points,
        shared_idx, shared_dist, running_dist,
    )
    # Note that this assumes first and last point are the same 
    if isinf(mid_x)
        mid_x, mid_y, running_dist = check_for_edge_mid(c1, 1, gap_idx - 1,
            shared_idx, shared_dist, running_dist,
        )
    end
    return mid_x, mid_y
end

"""
    points_in_poly(xy, coords::PolyVec{<:AbstractFloat})

Determines if the provided points are within the given polygon, including
checking that the points are not in any holes.
Inputs:
    xy  <Matrix{Real}> n-by-2 matrix of element where each row is a point and
        the first column is the x-coordinates and the second is y-coordinates
    coords  <PolyVec{AbstractFloat}> coordinates of polygon, with the exterior
        coordinates as the first element of the vector, any any hole coordinates
        as subsequent entries.
Outputs:
    in_idx  <Vector{Bool}> vector of booleans the length of the given points xy 
        where an entry is true if the corresponding element in xy is within the
        given polygon.
"""
function points_in_poly(xy, coords::PolyVec{<:AbstractFloat})
    in_idx = fill(false, length(xy[:, 1]))
    if !isempty(xy) && !isempty(coords[1][1])
        # Loop over exterior coords and each hole
        for i in eachindex(coords)
            in_on = inpoly2(xy, reduce(hcat, coords[i])')
            if i == 1  # Exterior outline of polygon - points must be within
                in_idx = in_idx .|| (in_on[:, 1] .|  in_on[:, 2])
            else  # Holes in polygon - points can't be within
                in_idx = in_idx .&& .!(in_on[:, 1] .|  in_on[:, 2])
            end
        end
    end
    return in_idx
end

"""
    points_in_poly(xy, multi_coords::Vector{<:PolyVec{<:AbstractFloat}})

Determines if the provided points are within the given multipolygon, including
checking that the points are not in any holes of any of the polygons.
Inputs:
    xy  <Matrix{Real}> n-by-2 matrix of element where each row is a point and
        the first column is the x-coordinates and the second is y-coordinates
    coords  <Vector{PolyVec{AbstractFloat}}> coordinates of the multi-polygon,
        with each element of the vector being a PolyVec of coordinates for a
        polygon and within each polygon the exterior coordinates as the first
        element of the PolyVec, any any hole coordinates as subsequent entries.
Outputs:
    in_idx  <Vector{Bool}> vector of booleans the length of the given points xy 
        where an entry is true if the corresponding element in xy is within the
        given polygon.
"""
function points_in_poly(xy, multi_coords::Vector{<:PolyVec{<:AbstractFloat}})
    # Check which of the points are within the domain coords
    in_idx = fill(false, length(xy[:, 1]))
    if !isempty(xy) && !isempty(multi_coords[1][1][1])
        # Loop over every polygon
        for i in eachindex(multi_coords)
            # See if the points are within current polygon
            in_idx = in_idx .|| points_in_poly(xy, multi_coords[i])
        end
    end
    return in_idx
end