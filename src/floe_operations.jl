"""
Translate each of the given coodinates by given vector -
Coordinates and vector must be vectors of same underlying type
Inputs: coords <Vector{Vector{Vector{Real}}}>
        vec <Vector{Real}>
Output: <Vector{Vector{Vector{Real}}}>
"""
function translate(coords::Vector{Vector{Vector{T}}}, vec::Vector{T}) where {T<:Real}
    return [c .+ repeat([vec], length(c)) for c in coords]
end

"""
Translate the given polygon by the given vector and return a new polygon -
Inputs: coords <LibGEOS.Polygon>
        vec <Vector{Real}>
Output: <LibGEOS.Polygon>
"""
function translate(poly::LG.Polygon, vec::Vector{T}) where {T<:Real}
    coords = LG.GeoInterface.coordinates(poly)::PolyVec64
    return LG.Polygon(translate(coords, convert(Vector{Float64}, vec)))
end

"""
Scale given polygon with respect to the reference point (0,0) - scaling factor is applied to both the x and y directions.
Inputs: coords <LibGEOS.Polygon>
        factor <Real>
Output: <LibGEOS.Polygon>
"""
function scale(poly::LG.Polygon, factor::T) where {T<:Real}
    coords = LG.GeoInterface.coordinates(poly)::PolyVec64
    return LG.Polygon(coords .* factor)
end


# Might be benefitial to add functions for coordinates
"""
Determine if polygon has one or more holes
Inputs: poly <LibGEOS.Polygon>
Outputs: <Bool>
"""
function hashole(poly::LG.Polygon)
    return LG.numInteriorRings(poly.ptr) > 0
end 

"""
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
Returns given polygon within a vector as it is the only region
Inputs: poly <LibGEOS.Polygon>
Outputs: <Vector{LibGEOS.Polygon}> 
"""
function sortregions(poly::LG.Polygon)
    return [poly]
end

"""
Sorts polygons within a multi-polygon by area in descending order
Inputs: multipoly <LibGEOS.MultiPolygon>
Outputs: <Vector{LibGEOS.Polygon}> 
"""
function sortregions(multipoly::LG.MultiPolygon)
    poly_lst = LG.getGeometries(multipoly)::Vector{LG.Polygon}
    return sort(poly_lst, by=LG.area, rev=true)
end

