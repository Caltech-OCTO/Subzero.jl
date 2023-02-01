
"""
    valid_ringvec(coords::RingVec{FT})

Takes a RingVec object and make sure that the last element has the same first element as
last element and that other than these two elements there are no duplicate, adjacent vertices.
Also asserts that the ring as at least three elements or else it cannot be made into a valid ring as it is a line segment. 
"""
function valid_ringvec!(ring::RingVec{FT}) where {FT<:AbstractFloat}
    deleteat!(ring, findall(i->ring[i]==ring[i+1], collect(1:length(ring)-1)))
    @assert length(ring) > 2 "Polgon needs at least 3 distinct points."
    if ring[1] != ring[end]
        push!(ring, ring[1])
    end
    return ring
end

"""
    valid_polyvec(coords::PolyVec{FT})

Takes a PolyVec object and make sure that the last element of each "ring" (vector of vector of floats)
has the same first element as last element and has not duplicate adjacent elements. Also asserts that each "ring"
as at least three distinct elements or else it is not a valid ring, but rather a line segment. 
"""
function valid_polyvec!(coords::PolyVec{FT}) where {FT<:AbstractFloat}
    for ring in coords
        valid_ringvec!(ring)
    end
    return coords
end

"""
    translate(coords::PolyVec{T}, vec

Translate each of the given coodinates by given vector -
Coordinates and vector must be vectors of same underlying type
Inputs: coords PolyVec{Float}
        vec <Vector{Real}>
Output: PolyVec{Float}
"""
function translate(coords::PolyVec{T}, vec) where {T<:AbstractFloat}
    return [c .+ repeat([vec], length(c)) for c in coords]
end

"""
    translate(poly::LG.Polygon, vec)

Translate the given polygon by the given vector and return a new polygon -
Inputs: coords <LibGEOS.Polygon>
        vec <Vector{Float}>
Output: <LibGEOS.Polygon>
"""
function translate(poly::LG.Polygon, vec)
    coords = LG.GeoInterface.coordinates(poly)::PolyVec{Float64}
    return LG.Polygon(translate(coords, convert(Vector{Float64}, vec)))
end

"""
    scale(poly::LG.Polygon, factor)

Scale given polygon with respect to the reference point (0,0).
Scaling factor is applied to both the x and y directions.
Inputs: coords <LibGEOS.Polygon>
        factor <Real>
Output: <LibGEOS.Polygon>
"""
function scale(poly::LG.Polygon, factor)
    coords = LG.GeoInterface.coordinates(poly)::PolyVec{Float64}
    return LG.Polygon(coords .* factor)
end

"""
    hashole(coords::PolyVec{FT})

Determine if polygon coordinates have one or more holes
Inputs: coords <PolyVec{Float}>
Outputs: <Bool>
"""
function hashole(coords::PolyVec{FT}) where FT<:AbstractFloat
    return length(coords) > 1
end

"""
    hashole(poly::LG.Polygon)

Determine if polygon has one or more holes
Inputs: poly <LibGEOS.Polygon>
Outputs: <Bool>
"""
function hashole(poly::LG.Polygon)
    return LG.numInteriorRings(poly.ptr) > 0
end 

"""
    hashole(multipoly::LG.MultiPolygon)

Determine if any of multipolygon's internal polygons has holes
Inputs: multipoly <LibGEOS.MultiPolygon>
Outputs: <Bool>
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
Inputs: coords <PolyVec{Float}>
Outputs: <PolyVec{Float}> 
"""
function rmholes(coords::PolyVec{FT}) where {FT<:AbstractFloat}
    return [coords[1]]
end

"""
    rmholes(poly::LG.Polygon)

Remove polygon's holes if they exist
Inputs: poly <LibGEOS.Polygon>
Outputs: <LibGEOS.Polygon> 
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
Inputs: multipoly <LibGEOS.MultiPolygon>
Outputs: <LibGEOS.MultiPolygon> 
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
Inputs: poly <LibGEOS.Polygon>
Outputs: <Vector{LibGEOS.Polygon}> 
"""
function sortregions(poly::LG.Polygon)
    return [poly]
end

"""
    sortregions(multipoly::LG.MultiPolygon)

Sorts polygons within a multi-polygon by area in descending order
Inputs: multipoly <LibGEOS.MultiPolygon>
Outputs: <Vector{LibGEOS.Polygon}> 
"""
function sortregions(multipoly::LG.MultiPolygon)
    poly_lst = LG.getGeometries(multipoly)::Vector{LG.Polygon}
    return sort(poly_lst, by=LG.area, rev=true)
end

"""
    seperate_xy(coords::PolyVec{T})

Pulls x and y coordinates from standard polygon vector coordinates into seperate vectors.
Only keeps external coordinates and disregards holes
Inputs: coords  PolyVec{Float} polygon coordinates
Outputs: x      <Vector{Float}> x coordinates
         y      <Vector{Float}> y coordinates
"""
function seperate_xy(coords::PolyVec{T}) where {T<:AbstractFloat}
    x = first.(coords[1])
    y = last.(coords[1])
    return x, y 
end

"""
    calc_moment_inertia(coords::PolyVec{T}, h; rhoice = 920.0)

Calculate the mass moment of intertia from polygon coordinates.
Inputs: coords      <PolyVec{Float}>
        h           <Real> height of floe
        rhoice      <Real> Density of ice
Output: mass moment of inertia <Float>
Note: Assumes that first and last point within the coordinates are the same. Will not give correct answer otherwise.

Based on paper: Marin, Joaquin."Computing columns, footings and gates through moments of area." Computers & Structures 18.2 (1984): 343-349.
"""
function calc_moment_inertia(coords::PolyVec{T}, h; ρi = 920.0) where {T<:AbstractFloat}
    x, y = seperate_xy(coords)
    N = length(x)
    wi = x[1:N-1] .* y[2:N] - x[2:N] .* y[1:N-1];
    Ixx = 1/12 * sum(wi .* ((y[1:N-1] + y[2:N]).^2 - y[1:N-1] .* y[2:N]))
    Iyy = 1/12 * sum(wi .* ((x[1:N-1] + x[2:N]).^2 - x[1:N-1] .* x[2:N]))
    return abs(Ixx + Iyy)*h*ρi;
end

"""
    calc_moment_inertia(poly::LG.Polygon, h; rhoice = 920.0)

Calculate the mass moment of intertia from a LibGEOS polygon object using above coordinate-based moment of intertia function.
Inputs: poly      LibGEOS.Polygon
        h           <Real> height of floe
        rhoice      <Real> Density of ice
Output: mass moment of inertia <Float>
"""
function calc_moment_inertia(poly::LG.Polygon, h; ρi = 920.0)
    return calc_moment_inertia(LG.GeoInterface.coordinates(poly),
                               h, ρi = 920.0)
end

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
Note: See note on calc_poly_angles for credit for this function.
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
    convex_angle_test(coords::RingVec{T}, t::Type{T} = Float64)

Determine which angles in the polygon are convex, with the assumption that the first angle is convex, no other
vertex has a smaller x-coordinate, and the vertices are assumed to be ordered in a clockwise sequence.
The test is based on the fact that every convex vertex is on the positive side of the line passing through
the two vertices immediately following each vertex being considered. 
Inputs:
        coords <RingVec{Float}> Vector of [x, y] vectors that make up the exterior of a polygon
        t  <AbstractFloat> datatype to run model with - must be a Float!
Outputs:
        sgn <Vector of 1s and -1s> One element for each [x,y] pair - if 1 then the angle at that vertex is convex,
        if it is -1 then the angle is concave.
"""
function convex_angle_test(coords::RingVec{T}, ::Type{T} = Float64) where T
    valid_ringvec!(coords)
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
        The positive side of the line should be in the right side of the vector (p2- p3).
        Δx and Δy give the direction of travel, establishing which of the extreme points (see above) should be on the + side.
        If that point is on the negative side of the line, then w is replaced by -w. =#
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
    calc_poly_angles(coords::PolyVec{T}, ::Type{T} = Float64))

Computes internal polygon angles (in degrees) of an arbitrary polygon given the coordinates ordered in a clockwise manner.
The program eliminates duplicate points, except that the first row must equal the last, so that the polygon is closed.
Inputs:
        coords  <PolyVec{Float}> coordinates from a polygon
        t       <AbstractFloat> datatype to run model with - must be a Float!
Outputs:
        Vector of polygon's interior angles

Note - Translated into Julia from the following program (including helper functions convex_angle_test and polyedge):
Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
Digital Image Processing Using MATLAB, Prentice-Hall, 2004
Revision: 1.6 Date: 2003/11/21 14:44:06
"""
function calc_poly_angles(coords::PolyVec{T}, ::Type{T} = Float64) where {T<:AbstractFloat}
    ext = valid_polyvec!(coords)[1]
    # Calculate needed vectors
    pdiff = diff(ext)
    npoints = length(pdiff)
    v1 = -pdiff
    v2 = vcat(pdiff[2:end], pdiff[1:1])
    v1_dot_v2 = [sum(v1[i] .* v2[i]) for i in collect(1:npoints)]
    mag_v1 = sqrt.([sum(v1[i].^2) for i in collect(1:npoints)])
    mag_v2 = sqrt.([sum(v2[i].^2) for i in collect(1:npoints)])
    # Protect against division by 0 caused by very close points
    mag_v1[mag_v1 .== 0.0] .= eps()
    mag_v2[abs.(mag_v2) .== 0.0] .= eps()
    angles = real.(acos.(v1_dot_v2 ./ mag_v1 ./ mag_v2) * 180 / pi)

    # The first angle computed was for the second vertex, and the last was for the first vertex.
    # Scroll one position down to make the last vertex be the first.
    sangles = circshift(angles, 1)
    # Now determine if any vertices are concave and adjust the angles accordingly.
    sgn = convex_angle_test(ext, T)
    sangles = [sgn[i] == -1 ? 360 - sangles[i] : sangles[i] for i in collect(1:length(sangles))]
    return sangles
end

"""
    calc_point_poly_dist(xp::Vector{T},yp::Vector{T}, vec_poly::PolyVec{T})

Compute the distances from each one of a set of np points on a 2D plane to a polygon.
Distance from point j to an edge k is defined as a distance from this point to a straight line passing
through vertices v(k) and v(k+1), when the projection of point j on this line falls INSIDE segment k;
and to the closest of v(k) or v(k+1) vertices, when the projection falls OUTSIDE segment k.
Inputs:
        xp  <Vector{Float}> x-coordinates of points to find distance from vec_poly
        yp  <Vector{Float}> y-coordiantes of points to find distance from vec_poly
        vec_poly    <PolyVec{Float}> coordinates of polygon
Outputs:
        List of distances from each point to the polygon. If the point is inside of the polygon the
        value will be negative.

Note - Translated into Julia from the following program:
p_poly_dist by Michael Yoshpe - last updated in 2006.
We mimic version 1 functionality with 4 inputs and 1 output.
Only needed code was translated.
"""
function calc_point_poly_dist(xp::Vector{T},yp::Vector{T}, vec_poly::PolyVec{T}) where {T<:AbstractFloat}
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
        min_dist, min_idx = findmin(p2c_dist, dims = 2)  # minimum distance to vertices
        # Coordinates of consecutive vertices
        V1 = Pv[1:end-1, :]
        V2 = Pv[2:end, :]
        Δv = V2 .-  V1
        # Vector of distances between each pair of consecutive vertices
        vds = hypot.(Δv[:, 1], Δv[:, 2])

        if (cumsum(vds)[end-1] - vds[end]) < 10eps(T)
            throw(ArgumentError("Polygon vertices should not lie on a straight line"))
        end

        # Each pair of consecutive vertices V1[j], V2[j] defines a rotated coordinate system
        # with origin at V1[j], and x axis along the vector V2[j]-V1[j].
        # Build the rotation matrix Cer from original to rotated system
        cθ = Δv[:, 1] ./ vds
        sθ = Δv[:, 2] ./  vds
        Cer = zeros(T, 2, 2, nv-1)
        Cer[1, 1, :] .= cθ
        Cer[1, 2, :] .= sθ
        Cer[2, 1, :] .= -sθ
        Cer[2, 2, :] .= cθ

        # Build the origin translation vector P1r in rotated frame by rotating the V1 vector
        V1r = hcat(cθ .* V1[:, 1] .+ sθ .* V1[:, 2],
                -sθ .* V1[:, 1] .+ cθ .* V1[:, 2])

        # Ppr is a 3D array of size 2*np*(nv-1). Ppr(1,j,k) is an X coordinate of point
        # j in coordinate systems defined by segment k. Ppr(2,j,k) is its Y coordinate.
        Ppr = zeros(T, 2, np, nv-1)
        # Rotation and Translation
        Ppr[1, :, :] .= Pp * Cer[1, :, :] .-
                        permutedims(repeat(V1r[:, 1], 1, 1, np), [2, 3, 1])[1, :, :]
        Ppr[2, :, :] .= Pp * Cer[2, :, :] .-
                        permutedims(repeat(V1r[:, 2], 1, 1, np), [2, 3, 1])[1, :, :]

        # x and y coordinates of the projected (cross-over) points in original
        # coordinate
        r = Ppr[1, :, :]
        cr = Ppr[2, :, :]
        B = fill(convert(T, Inf), np, nv-1)
        # For the projections that fall inside the segments, find the minimum distances from points to
        # their projections (note, that for some points these might not exist)
        for i in eachindex(r)
            if r[i] > 0 && r[i] < vds[cld(i, np)]
                B[i] = cr[i]
            end
        end
        cr_min, cr_min_idx = findmin(abs.(B), dims = 2)
        # For projections that fall outside segments, closest point is a vertex
        # These points have a negative value if point is actually outside of polygon
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
        T[]
    end
    return min_lst
end

"""
    intersect_lines(l1, l2)

Finds the intersection points of two curves l1 and l2. The curves l1, l2 can be either closed or open. In this version, l1 and l2 must be distinct. If no intersections are found, the returned P is empty.
Inputs:
        l1 <PolyVec{Float}> line/polygon coordinates
        l2 <PolyVec{Float}> line/polygon coordinates
Outputs:
        N intersection points in a Nx2 matrix where column 1 is the x-coordinates and column 2 is the y-coordinates and each intersection point is a row.

Note - Translated into Julia from the following program:
NS (2022). Curve intersections (https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections), MATLAB Central File Exchange. Retrieved November 2, 2022.
Only translated for the case where l1 and l2 are distinct. 
"""
function intersect_lines(l1, l2)
    x1, y1 = seperate_xy(l1)
    x2, y2 = seperate_xy(l2)
    x2t = x2'
    y2t = y2'
    Δx1 = diff(x1)
    Δx2 = diff(x2t, dims = 2)
    Δy1 = diff(y1)
    Δy2 = diff(y2t, dims = 2)

    # Determine 'signed distances' 
    S1 = Δx1 .* y1[1:end-1] - Δy1 .* x1[1:end-1]
    s1 = (Δx1 .* y2t .- Δy1 .*x2t)  # Needed for S1 calculation
    C1 = (s1[:, 1:end-1] .- S1) .* (s1[:, 2:end] .- S1) .<= 0

    S2 = Δx2 .* y2t[:, 1:end-1] - Δy2 .* x2t[:, 1:end-1]
    s2 = (y1 .* Δx2 .- x1 .* Δy2)'  # Needed for S2 calculation
    C2 = ((s2[:, 1:end-1] .- S2') .* (s2[:, 2:end] .- S2') .<= 0)'

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

function cutpolygon_coords(poly_coords, yp, ::Type{T} = Float64) where T
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
    # Add unique points to coordinate list for new polygon
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
    # If no NaNs, just add to list
    if isempty(nanidx_all)
        if new_poly_coords[1] != new_poly_coords[end]
            push!(new_poly_coords, new_poly_coords[1])
        end
        if length(new_poly_coords) > 3
            push!(new_polygons, [new_poly_coords])
        end
    # Seperate out NaNs to seperate out polygons multiple polygons and add to list
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
Assumes that polygon provided has a hole! If not, it will error.
"""
function split_polygon_hole(poly::LG.Polygon, ::Type{T} = Float64) where T
    poly_coords = LG.GeoInterface.coordinates(poly)
    full_coords = [poly_coords[1]]  # Without any holes
    h1 = LG.Polygon([poly_coords[2]])
    h1_center = LG.GeoInterface.coordinates(LG.centroid(h1))
    poly_bottom = LG.MultiPolygon(cutpolygon_coords(full_coords, h1_center[2], T))
    poly_bottom =  LG.intersection(poly_bottom, poly)
    poly_top = LG.difference(poly, poly_bottom)
    return LG.getGeometries(poly_bottom), LG.getGeometries(poly_top)
end