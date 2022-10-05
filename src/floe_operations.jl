
"""
    valid_ringvec(coords::RingVec{FT})

Takes a RingVec object and make sure that the last element has the same first element as last element. Asserts that the ring as at least three elements or else it cannot be made into a valid ring as it is a line segment. 
"""
function valid_ringvec!(ring::RingVec{FT}) where {FT<:AbstractFloat}
    if ring[1] != ring[end]
        @assert length(ring) > 2 "Polgon needs at least 3 points."
        push!(ring, ring[1])
    end
    return ring
end

"""
    valid_polyvec(coords::PolyVec{FT})

Takes a PolyVec object and make sure that the last element of each "ring" (vector of vector of floats) has the same first element as last element. Asserts that each "ring" as at least three elements or else it is not a valid ring, but rather a line segment. 
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

Scale given polygon with respect to the reference point (0,0) - scaling factor is applied to both the x and y directions.
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

Pulls x and y coordinates from standard polygon vector coordinates into seperate vectors - only keeps external coordinates and disregards holes
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