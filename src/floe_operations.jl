function hashole(poly::LibGEOS.Polygon)
    return LibGEOS.numInteriorRings(poly.ptr) > 0
end 

function hashole(multipoly::LibGEOS.MultiPolygon)
    poly_lst = LibGEOS.getGeometries(multipoly)
    for poly in poly_lst
        if hashole(poly)
            return true
        end
    end
    return false
end