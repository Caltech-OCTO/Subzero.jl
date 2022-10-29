
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
function convex_angle_test(coords::RingVec{T}, t::Type{T} = Float64) where T
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
    calc_poly_angles(coords::PolyVec{T})

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
function calc_poly_angles(coords::PolyVec{T}, t::Type{T} = Float64) where {T<:AbstractFloat}
    ext = valid_polyvec!(coords)[1]
    # Calculate needed vectors
    pdiff = ext[2:end] .- ext[1:end-1]
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
    circshift!(angles, -1)
    # Now determine if any vertices are concave and adjust the angles accordingly.
    sgn = convex_angle_test(ext, T)
    angles = [sgn[i] == -1 ? 360 - angles[i] : angles[i] for i in collect(1:length(angles))]
    return angles
end

function calc_poly_dist(coords::PolyVec{T}) where {T<:AbstractFloat}
    return
end