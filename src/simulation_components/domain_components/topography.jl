export TopographyElement  # Topographic element within domain

"""
    TopographyElement{FT}<:AbstractDomainElement{FT}

Singular topographic element with coordinates field storing where the element is
within the grid. These are used to create the desired topography within the
simulation and will be treated as islands or coastline within the model
in that they will not move or break due to floe interactions, but they will
affect floes.
"""
struct TopographyElement{FT}<:AbstractDomainElement{FT}
    poly::Polys{FT}
    centroid::Vector{FT}
    rmax::FT

    function TopographyElement{FT}(
        poly,
        centroid,
        rmax,
    ) where {FT <: AbstractFloat}
        if rmax <= 0
            throw(ArgumentError("Topography element maximum radius must be \
                positive and non-zero."))
        end
        poly = GO.ClosedRing()(poly)
        new{FT}(poly, centroid, rmax)
    end
end

function TopographyElement(::Type{FT} = Float64; poly = nothing, coords = nothing) where FT
    if isnothing(poly) && !isnothing(coords)
        poly = make_polygon(coords)
    elseif isnothing(poly)
        throw(ArgumentError("To create a topography element the user must provide either a polygon (with the poly keyword) or coordinates (with the coord keyword)."))
    end
    # Clean up polygon and calculate centroid and maximum radius
    poly = GO.ClosedRing()(poly)
    rmholes!(poly)
    centroid = collect(GO.centroid(poly)) # TODO: Remove collect once type is changed
    rmax = calc_max_radius(poly, centroid, FT)
    return TopographyElement{FT}(poly, centroid, rmax)
end

"""
    initialize_topography_field(
        ::Type{FT},
        coords,
    )

Create a field of topography from a list of polygon coordiantes.
Inputs:
    Type{FT}        <AbstractFloat> Type for grid's numberical fields -
                        determines simulation run type
    coords          <Vector{PolyVec}> list of polygon coords to make into floes
Outputs:
    topo_arr <StructArray{TopographyElement}> list of topography elements
    created from given polygon coordinates
"""
function initialize_topography_field(
    ::Type{FT} = Float64; coords
) where {FT <: AbstractFloat}
    topo_multipoly = GO.DiffIntersectingPolygons()(GI.MultiPolygon(coords))
    topo_arr = StructArray{TopographyElement{FT}}(undef, GI.npolygon(topo_multipoly))
    for (i, poly) in enumerate(GI.getpolygon(topo_multipoly))
        topo_arr[i] = TopographyElement(FT; poly)
    end
    return topo_arr
end