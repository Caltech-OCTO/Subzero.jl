
"""
    valid_ringvec(coords::RingVec{FT})

Takes a RingVec object and make sure that the last element has the same first
element as last element and that other than these two elements there are no
duplicate, adjacent vertices. Also asserts that the ring as at least three
elements or else it cannot be made into a valid ring as it is a line segment. 
"""
function valid_ringvec!(ring::RingVec{FT}) where {FT<:AbstractFloat}
    deleteat!(ring, findall(i->ring[i]==ring[i+1], collect(1:length(ring)-1)))
    if ring[1] != ring[end]
        push!(ring, ring[1])
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
function valid_polyvec!(coords::PolyVec{FT}) where {FT<:AbstractFloat}
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
find_poly_coords(poly::LG.Polygon) =
    LG.GeoInterface.coordinates(poly)::PolyVec{Float64}
# [LG.GeoInterface.coordinates(LG.exteriorRing(poly))]::PolyVec{Float64}


"""
    find_multipoly_coords(poly)

Syntactic sugar for using LibGEOS to find a multipolygon's coordinates
Input:
    poly    <LibGEOS.Polygon>
Output:
    <Vector{PolyVec}> representing the floe's coordinates xy plane
"""
find_multipoly_coords(poly::LG.Polygon) =
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
    translate(coords::PolyVec{T}, vec

Translate each of the given coodinates by given vector -
Coordinates and vector must be vectors of same underlying type
Inputs:
    coords PolyVec{Float}
    vec <Vector{Real}>
Output:
    PolyVec{Float}
"""
function translate(coords::PolyVec{T}, vec) where {T<:AbstractFloat}
    return [c .+ repeat([vec], length(c)) for c in coords]
end

"""
    translate(poly::LG.Polygon, vec)

Translate the given polygon by the given vector and return a new polygon -
Inputs:
    coords <LibGEOS.Polygon>
    vec <Vector{Float}>
Output:
    <LibGEOS.Polygon>
"""
function translate(poly::LG.Polygon, vec)
    coords = find_poly_coords(poly)
    return LG.Polygon(translate(coords, convert(Vector{Float64}, vec)))
end

"""
    scale(poly::LG.Polygon, factor)

Scale given polygon with respect to the reference point (0,0).
Scaling factor is applied to both the x and y directions.
Inputs:
    coords <LibGEOS.Polygon> 
    factor <Real>
Output: <LibGEOS.Polygon>
"""
function scale(poly::LG.Polygon, factor)
    coords = find_poly_coords(poly)
    return LG.Polygon(coords .* factor)
end

"""
    rotate_radians(coords::PolyVec, α)

Rotate a polygon's coordinates by α radians around the origin.
Inputs:
    coords  <PolyVec{AbstractFloat}> polygon coordinates
    α       <Real> radians to rotate the coordinates
Outputs:
    <PolyVec{AbstractFloat}> coords rotates by α radians
"""
function rotate_radians(coords::PolyVec, α)
    return [map(p -> [cos(α)*p[1] - sin(α)*p[2],
                         sin(α)*p[1] + cos(α)p[2]], coords[1])]
end

"""
    rotate_degrees(coords::PolyVec, α)

Rotate a polygon's coordinates by α degrees around the origin.
Inputs:
    coords <PolyVec{AbstractFloat}> polygon coordinates
    α       <Real> degrees to rotate the coordinates
Outputs:
    <PolyVec{AbstractFloat}> coords rotates by α degrees
"""
rotate_degrees(coords::PolyVec, α) = rotate_radians(coords, α * π/180)

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
function hashole(poly::LG.Polygon)
    return LG.numInteriorRings(poly.ptr) > 0
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

"""
    rmholes(poly::LG.Polygon)

Remove polygon's holes if they exist
Inputs:
    poly <LibGEOS.Polygon> LibGEOS polygon
Outputs:
    <LibGEOS.Polygon>  LibGEOS polygon without any holes
"""
function rmholes(poly::LG.Polygon)
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
function sortregions(poly::LG.Polygon)
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
    poly_lst = LG.getGeometries(multipoly)::Vector{LG.Polygon}
    return sort(poly_lst, by=LG.area, rev=true)
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
calc_moment_inertia(poly::LG.Polygon, h; ρi = 920.0) = 
    calc_moment_inertia(
        find_poly_coords(poly),
        find_poly_centroid(poly),
        h,
        ρi = ρi,
    )

"""
    polyedge(p1, p2, t::Type{T} = Float64)

Outputs the coefficients of the line passing through p1 and p2.
The line is of the form w1x + w2y + w3 = 0. 
Inputs:
    p1 <Vector{Float}> [x, y] point
    p2 <Vector{Float}> [x, y] point
    t  <AbstractFloat> datatype to run model with - must be a Float!
Outputs:
    Three-element vector for coefficents of line passing through p1 and p2
Note:
    See note on calc_poly_angles for credit for this function.
"""
function polyedge(p1, p2, t::Type{T} = Float64) where T 
    x1 = p1[1]
    y1 = p1[2]
    x2 = p2[1]
    y2 = p2[2]
    w = if x1 == x2
            [-1/x1, 0, 1]
        elseif y1 == y2
            [0, -1/y1, 1]
        elseif x1 == y1 && x2 == y2
            [1, 1, 0]
        else
            v = (y1 - y2)/(x1*(y2 - y1) - y1*(x2 - x1) + eps(T))
            [v, -v*(x2 - x1)/(y2 - y1), 1]
        end
    return w
end

"""
    orient_coords(coords::RingVec{T}) where T

Take given coordinates and make it so that the first point has the smallest
x-coordiante and so that the coordinates are ordered in a clockwise sequence.
Duplicates vertices will be removed and the coordiantes will be closed (first
and last point are the same).

Input:
    coords  <RingVec> vector of points [x, y]
Output:
    coords  <RingVec> oriented clockwise with smallest x-coordinate first
"""
function orient_coords(coords::RingVec{T}) where T
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
    coords = circshift(coords, -extreem_idx + 1)
    valid_ringvec!(coords)

    # if coords are counterclockwise, switch to clockwise
    orient_matrix = hcat(
        ones(T, 3),
        vcat(coords[1]', coords[2]', coords[end-1]') # extreem and adjacent points
    )
    if det(orient_matrix) > 0
        reverse!(coords)
    end
    return coords
end

"""
    convex_angle_test(coords::RingVec{T})

Determine which angles in the polygon are convex, with the assumption that the
first angle is convex, no other vertex has a smaller x-coordinate, and the
vertices are assumed to be ordered in a clockwise sequence. The test is based on
the fact that every convex vertex is on the positive side of the line passing
through the two vertices immediately following each vertex being considered. 
Inputs:
    coords <RingVec{Float}> Vector of [x, y] vectors that make up the exterior
        of a polygon
Outputs:
        sgn <Vector of 1s and -1s> One element for each [x,y] pair - if 1 then
            the angle at that vertex is convex, if it is -1 then the angle is
            concave.
"""
function convex_angle_test(coords::RingVec{T}) where T
    L = 10^25
    # Extreme points used in following loop, apended by a 1 for dot product
    top_left = [-L, -L, 1]
    top_right = [-L, L, 1]
    bottom_left = [L, -L, 1]
    bottom_right = [L, L, 1]
    sgn = [1]  # First vertex is convex

    for k in collect(2:length(coords)-1)
        p1 = coords[k-1]
        p2 = coords[k]  # Testing this point for concavity
        p3 = coords[k+1]
        # Coefficents of polygon edge passing through p1 and p2
        w = polyedge(p1, p2)

        #= Establish the positive side of the line w1x + w2y + w3 = 0.
        The positive side of the line should be in the right side of the vector
        (p2- p3).Δx and Δy give the direction of travel, establishing which of
        the extreme points (see above) should be on the + side. If that point is
        on the negative side of the line, then w is replaced by -w. =#
        Δx = p2[1] - p1[1]
        Δy = p2[2] - p1[2]
        if Δx == Δy == 0
            throw(ArgumentError("Data into convextiy test is 0 or duplicated"))
        end
        vector_product =
            if Δx <= 0 && Δy >= 0  # Bottom_right should be on + side.
                dot(w, bottom_right)
            elseif Δx <= 0 && Δy <=0  # Top_right should be on + side.
                dot(w, top_right)
            elseif Δx>=0 && Δy<=0  # Top_left should be on + side.
                dot(w, top_left)
            else  # Bottom_left should be on + side.
                dot(w, bottom_left)
            end
            w *= sign(vector_product)

            # For vertex at p2 to be convex, p3 has to be on + side of line
            if (w[1]*p3[1] + w[2]*p3[2] + w[3]) < 0
                push!(sgn, -1)
            else
                push!(sgn, 1)
            end
    end
    return sgn
end

"""
    calc_poly_angles(coords::PolyVec{T})

Computes internal polygon angles (in degrees) of an arbitrary simple polygon.
The program eliminates duplicate points, except that the first row must equal
the last, so that the polygon is closed.
Inputs:
    coords  <PolyVec{Float}> coordinates from a polygon
Outputs:
    Vector of polygon's interior angles in degrees

Note - Translated into Julia from the following program (including helper
    functions convex_angle_test and polyedge):
    Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
    Digital Image Processing Using MATLAB, Prentice-Hall, 2004
    Revision: 1.6 Date: 2003/11/21 14:44:06
"""
function calc_poly_angles(coords::PolyVec{T}) where {T<:AbstractFloat}
    ext = orient_coords(coords[1]) # ignore any holes in the polygon
    # Calculate needed vectors
    pdiff = diff(ext)
    npoints = length(pdiff)
    v1 = -pdiff
    v2 = vcat(pdiff[2:end], pdiff[1:1])
    v1_dot_v2 = [sum(v1[i] .* v2[i]) for i in collect(1:npoints)]
    mag_v1 = sqrt.([sum(v1[i].^2) for i in collect(1:npoints)])
    mag_v2 = sqrt.([sum(v2[i].^2) for i in collect(1:npoints)])
    # Protect against division by 0 caused by very close points
    replace!(mag_v1, 0=>eps(T))
    replace!(mag_v2, 0=>eps(T))
    angles = real.(acos.(v1_dot_v2 ./ mag_v1 ./ mag_v2) * 180 / pi)

    #= The first angle computed was for the second vertex, and the last was for
    the first vertex. Scroll one position down to make the last vertex be the
    first. =#
    sangles = circshift(angles, 1)
    # Now determine if any vertices are concave and adjust angles accordingly.
    sgn = convex_angle_test(ext)
    for i in eachindex(sangles)
        sangles[i] = (sgn[i] == -1) ? (-sangles[i] + 360) : sangles[i]
    end
    return sangles
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
    <Matrix{AbstractFloat}> N intersection points in a Nx2 matrix where column 1
        is the x-coordinates and column 2 is the y-coordinates and each row is
        an intersection point.

Note - Translated into Julia from the following program:
    NS (2022). Curve intersections
    (https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections),
    MATLAB Central File Exchange. Retrieved November 2, 2022.
    Only translated for the case where l1 and l2 are distinct. 
"""
function intersect_lines(l1, l2)
    x1, y1 = separate_xy(l1)
    x2, y2 = separate_xy(l2)
    x2t = x2'
    y2t = y2'
    Δx1 = diff(x1)
    Δx2 = diff(x2t, dims = 2)
    Δy1 = diff(y1)
    Δy2 = diff(y2t, dims = 2)

    # Determine 'signed distances' 
    S1 = Δx1 .* @view(y1[1:end-1]) .- Δy1 .* @view(x1[1:end-1])
    s1 = (Δx1 .* y2t .- Δy1 .*x2t)  # Needed for S1 calculation
    C1 = (@view(s1[:, 1:end-1]) .- S1) .* (@view(s1[:, 2:end]) .- S1) .<= 0

    S2 = Δx2 .* @view(y2t[:, 1:end-1]) .- Δy2 .* @view(x2t[:, 1:end-1])
    s2 = (y1 .* Δx2 .- x1 .* Δy2)'  # Needed for S2 calculation
    C2 = ((@view(s2[:, 1:end-1]) .- S2') .* (@view(s2[:, 2:end]) .- S2') .<= 0)'

    # Obtain the segments where an intersection is expected
    idx = findall(C1 .& C2)

    P = if isempty(idx)
        zeros(eltype(x1), 2,0)
    else
        # Transpose and prepare for output
        i = getindex.(idx, 1)
        j = getindex.(idx, 2)[:, :]
        Δx2t = Δx2'
        Δy2t = Δy2'
        S2t = S2'
        L = Δy2t[j] .* Δx1[i] - Δy1[i] .* Δx2t[j]
        i = i[:, :][L .!= 0]
        j = j[L .!= 0]
        L = L[L .!= 0]
        # Solve system of eqs to get the common points
        unique(hcat((Δx2t[j] .* S1[i] - Δx1[i] .* S2t[j]) ./ L,
                    (Δy2t[j] .* S1[i] - Δy1[i] .* S2t[j]) ./ L), dims = 1)
    end
    return P
end

"""
    cut_polygon_coords(poly_coords::PolyVec, yp, ::Type{T} = Float64)

Cut polygon through the line y = yp and return the polygon(s) coordinates below
the line
Inputs:
    poly_coords <PolyVec>   polygon coordinates
    yp          <Float>     value of line to split polygon through using line y = yp
                <Type>      Type of abstract float to run simulation with
Outputs:
    new_polygons <Vector{PolyVec}> List of coordinates of polygons below line y = yp. 
Note: Code translated from MATLAB to Julia. Credit for initial code to Dominik
    Brands (2010) and Jasper Menger (2009). Only needed pieces of function are
    translated (horizonal cut).
"""
function cut_polygon_coords(poly_coords::PolyVec, yp, ::Type{T} = Float64) where T
    # Loop through each edge
    coord1 = poly_coords[1][1:end-1]
    coord2 = poly_coords[1][2:end]
    for i in eachindex(coord1)
        x1, y1 = coord1[i]
        x2, y2 = coord2[i]
        # If both edge endpoints are above cut line, remove edge
        if y1 > yp && y2 > yp
            coord1[i] = [NaN, NaN]
            coord2[i] = [NaN, NaN]
        # If start point is above cut line, move down to intersection point
        elseif y1 > yp
            coord1[i] = [(yp - y2)/(y1 - y2) * (x1 - x2) + x2, yp]
        # If end point is above cut line, move down to intersection point
        elseif y2 > yp
            coord2[i] = [(yp - y1)/(y2 - y1) * (x2 - x1) + x1, yp]
        end
    end
    # Add non-repeat points to coordinate list for new polygon
    new_poly_coords = [coord1[1]]
    for i in eachindex(coord1)
        if !isequal(coord1[i], new_poly_coords[end])
            push!(new_poly_coords, coord1[i])
        end
        if !isequal(coord2[i], new_poly_coords[end])
            push!(new_poly_coords, coord2[i])
        end
    end

    new_polygons = Vector{PolyVec{T}}()
    # Multiple NaN's indicate new polygon if they seperate coordinates
    nanidx_all = findall(c -> isnan(sum(c)), new_poly_coords)
    # If no NaNs, just add coordiantes to list
    if isempty(nanidx_all)
        if new_poly_coords[1] != new_poly_coords[end]
            push!(new_poly_coords, new_poly_coords[1])
        end
        if length(new_poly_coords) > 3
            push!(new_polygons, [new_poly_coords])
        end
    # Seperate out NaNs to seperate multiple polygons and add to list
    else
        if new_poly_coords[1] == new_poly_coords[end]
            new_poly_coords = new_poly_coords[1:end-1]
        end
        # Shift so each polygon's vertices are together in a section
        new_poly_coords = circshift(new_poly_coords, -nanidx_all[1] + 1)
        nanidx_all .-= (nanidx_all[1] - 1)
        # Determine start and stop point for each polygon's coordinates
        start_poly = nanidx_all .+ 1
        end_poly = [nanidx_all[2:end] .- 1; length(new_poly_coords)]
        for i in eachindex(start_poly)
            if start_poly[i] <= end_poly[i]
                poly = new_poly_coords[start_poly[i]:end_poly[i]]
                if poly[1] != poly[end]
                    push!(poly, poly[1])
                end
                if length(poly) > 3
                    push!(new_polygons, [poly])
                end
            end
        end
    end

    return new_polygons
end

"""
    split_polygon_hole(poly::LG.Polygon, ::Type{T} = Float64)

Splits polygon horizontally through first hole and return lists of polygons
created by split.
Inputs:
        poly    <LG.Polygon> polygon to split
                <Type> Float type to run simulation with
Outputs:
    <(Vector{LibGEOS.Polyon}, (Vector{LibGEOS.Polyon}> list of polygons created
    from split through first hole below line and polygons through first hole
    above line. Note that if there is no hole, a list of the original polygon
    and an empty list will be returned
"""
function split_polygon_hole(poly::LG.Polygon, ::Type{T} = Float64) where T
    bottom_list = Vector{LG.Polygon}()
    top_list = Vector{LG.Polygon}()
    if hashole(poly)  # Polygon has a hole
        poly_coords = find_poly_coords(poly)
        full_coords = [poly_coords[1]]
        h1 = LG.Polygon([poly_coords[2]])  # First hole
        h1_center = find_poly_centroid(h1)
        poly_bottom = LG.MultiPolygon(
            cut_polygon_coords(full_coords, h1_center[2], T)
        )
         # Adds in any other holes in poly
        poly_bottom =  LG.intersection(poly_bottom, poly)
        poly_top = LG.difference(poly, poly_bottom)
        bottom_list = LG.getGeometries(poly_bottom)::Vector{LG.Polygon}
        top_list = LG.getGeometries(poly_top)::Vector{LG.Polygon}
    else  # No hole
        bottom_list, top_list = Vector{LG.Polygon}([poly]), Vector{LG.Polygon}()
    end
    return bottom_list, top_list
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