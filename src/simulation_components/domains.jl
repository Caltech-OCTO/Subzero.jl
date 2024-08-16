# Cardinal directions
export AbstractDirection, North, South, East, West
# All elements that make up the domain
export AbstractDomainElement
# Boundary of the model - property of each of the 4 walls (north, south, east, west)
export AbstractBoundary, OpenBoundary, PeriodicBoundary, CollisionBoundary, MovingBoundary
# Topographic element within domain
export TopographyElement
# Domain definition (combines 4 boundaries and topography)
export Domain

"""
    YourDirection <: AbstractDirection

Each domain within a Subzero.jl [`Model`](@ref) must have four (4) boundaries (subtypes
of [`AbstractBoundary`](@ref)) where each of these boundaries is parametrically typed by the
direction of the boundary. The user will first choose one of the four cardinal directions,
the subtypes of `AbstractDirection`:
- [`North`](@ref)
- [`South`](@ref)
- [`East`](@ref)
- [`West`](@ref)

This abstract type is not meant to be extended by the user, unless the user wants to move
away from a rectangular domain with assigned cardinal directions for each wall. This would
a more major redesign and the user should check out the [developer documentation]("devdocs.md").
"""
abstract type AbstractDirection end

"""
    North<:AbstractDirection

A simple subtype of [`AbstractDirection`](@ref) used for parametrically typing a subtype of
[`AbstractBoundary`](@ref) if that boundary is the northern boundary in a rectangular domain.
"""
struct North<:AbstractDirection end

#=
    _grid_boundary_info(FT, North, grid)

Create a bounding box polygon representing the Northern boundary of a rectangular grid with
points of float type `FT`. If the length of the grid x-extent is `Lx = xf - x0` then the
x-extent of the boundary polygon will range from `x0 - Lx/2` to `xf + Lx/2` in the x-direction.
If the length of the grid y-extent is `Ly = yf - y0` then the boundary polygon will range from
`yf` to `yf + Ly/2` in the y-direction. This will create overlap with other boundary walls,
if all created from the grid, making sure all floes connect with boundaries at edges of
the domain.

Also return value `yf` as the value representing the edge of the boundary connecting with
the edge of the domain.
=#
function _grid_boundary_info(::Type{FT}, ::Type{North}, grid::RegRectilinearGrid) where FT
    Δx = (grid.xf - grid.x0)/2
    Δy = (grid.yf - grid.y0)/2
    poly =  _make_bounding_box_polygon(FT, grid.x0 - Δx, grid.xf + Δx, grid.yf, grid.yf + Δy)
    return poly, grid.yf
end

"""
    South<:AbstractDirection

A simple subtype of [`AbstractDirection`](@ref) used for parametrically typing a subtype of
[`AbstractBoundary`](@ref) if that boundary is the southern boundary in a rectangular domain.
"""
struct South<:AbstractDirection end

#=
    _grid_boundary_info(FT, South, grid)

Create a bounding box polygon representing the Southern boundary of a rectangular grid with
points of float type `FT`. If the length of the grid x-extent is `Lx = xf - x0` then the
x-extent of the boundary polygon will range from `x0 - Lx/2` to `xf + Lx/2` in the x-direction.
If the length of the grid y-extent is `Ly = yf - y0` then the boundary polygon will range from
`y0 - Ly/2` to `y0` in the y-direction. This will create overlap with other boundary walls,
if all created from the grid, making sure all floes connect with boundaries at edges of
the domain.

Also return value `y0` as the value representing the edge of the boundary connecting with
the edge of the domain.
=#
function _grid_boundary_info(::Type{FT}, ::Type{South}, grid::RegRectilinearGrid) where FT
    Δx = (grid.xf - grid.x0)/2
    Δy = (grid.yf - grid.y0)/2
    poly = _make_bounding_box_polygon(FT, grid.x0 - Δx, grid.xf + Δx, grid.y0 - Δy, grid.y0)
    return poly, grid.y0
end

"""
    East<:AbstractDirection


A simple subtype of [`AbstractDirection`](@ref) used for parametrically typing a subtype of
[`AbstractBoundary`](@ref) if that boundary is the eastern boundary in a rectangular domain.
"""
struct East<:AbstractDirection end

#=
    _grid_boundary_info(FT, East, grid)

Create a bounding box polygon representing the Eastern boundary of a rectangular grid with
points of float type `FT`. If the length of the grid x-extent is `Lx = xf - x0` then the
x-extent of the boundary polygon will range from `xf` to `xf + Lx/2` in the x-direction.
If the length of the grid y-extent is `Ly = yf - y0` then the boundary polygon will range from
`y0 - Ly/2` to `yf + Ly/2` in the y-direction. This will create overlap with other boundary walls,
if all created from the grid, making sure all floes connect with boundaries at edges of
the domain.

Also return value `xf` as the value representing the edge of the boundary connecting with
the edge of the domain.
=#
function _grid_boundary_info(::Type{FT}, ::Type{East}, grid::RegRectilinearGrid) where FT
    Δx = (grid.xf - grid.x0)/2
    Δy = (grid.yf - grid.y0)/2
    poly = _make_bounding_box_polygon(FT, grid.xf, grid.xf + Δx, grid.y0 - Δy, grid.yf + Δy)
    return poly, grid.xf
end

"""
    West<:AbstractDirection


A simple subtype of [`AbstractDirection`](@ref) used for parametrically typing a subtype of
[`AbstractBoundary`](@ref) if that boundary is the western boundary in a rectangular domain.
"""
struct West<:AbstractDirection end

#=
    _grid_boundary_info(FT, West, grid)

Create a bounding box polygon representing the Western boundary of a rectangular grid with
points of float type `FT`. If the length of the grid x-extent is `Lx = xf - x0` then the
x-extent of the boundary polygon will range from `x0 - Lx/2` to `x0` in the x-direction.
If the length of the grid y-extent is `Ly = yf - y0` then the boundary polygon will range from
`y0 - Ly/2` to `yf + Ly/2` in the y-direction. This will create overlap with other boundary walls,
if all created from the grid, making sure all floes connect with boundaries at edges of
the domain.

Also return value `x0` as the value representing the edge of the boundary connecting with
the edge of the domain.
=#
function _grid_boundary_info(::Type{FT}, ::Type{West}, grid::RegRectilinearGrid) where FT
    Δx = (grid.xf - grid.x0)/2
    Δy = (grid.yf - grid.y0)/2
    poly = _make_bounding_box_polygon(FT, grid.x0 - Δx, grid.x0, grid.y0 - Δy, grid.yf + Δy)
    return poly, grid.x0
end

"""
    YourDomainElement{FT} <: AbstractDomainElement{FT}

An abstract type for all of the element that create the shape and behavior of the [`Domain`](@ref):
the 4 boundary walls that make up the rectangular domain and the topography within the `Domain`. 

When building a [`Domain`](@ref), the user will create four (4) boundary elements that are
concrete types of [`AbstractBoundary`](@ref), which is itself a subtype of
`AbstractDomainElement`, and a list of [`TopographyElement`](@ref) objects, which are also a
concrete subtype of `AbstractDomainElement`.

## _API_
If the user wanted to add more domain elements, other than boundaries and topography, they
would create a new subtype of `AbstractDomainElement` and would need to write methods for
the following functions:
- 

If the user then wanted these in the [`Domain`](@ref), the user would need to change the
`Domain` struct.
"""
abstract type AbstractDomainElement{FT<:AbstractFloat} end

"""
    AbstractBoundary{D, FT} <: AbstractDomainElement{FT}

When running a Subzero simulation, each simulation has a rectangular domain that contains
the ice pack, made up for four boundaries, one in each of the cardinal directions:
[`North`](@ref), [`South`](@ref), [`East`](@ref), and [`West`](@ref). The user will pick
four boundary types, each concrete subtypes of `AbstractBoundary` and assign each of them a
direction, `D` and a float type `FT`, by typing them [parametrically]("https://docs.julialang.org/en/v1/manual/types/#Parametric-Types"). 

The user can pick from the currently implemented concrete subtypes of `AbstractBoundary`
([`OpenBoundary`](@ref), [`CollisionBoundary`](@ref), [`PeriodicBoundary`](@ref), and
[`MovingBoundary`](@ref)) or implement their own.

For now, each boundary wall is assumed to be a rectangle. For correct behavior, the corners
of each of the boundary rectangles should overlap to make sure that no floes reach outside
of the domain without touching a boundary. These rectangles should be represented as
polygons within a `poly` field within each boundary. Furthermore, each boundary should have
a field `val` that represents either the equation `y = val` (for `North` and `South` walls)
or `x = val` (for `East` and `West` walls) that denotes the line that marks the "inner-most"
edge of the domain. If an ice floe passes over the boudnary's line, the floe interacts with
the boundary.

```text
 ________________
|__|____val___|__| <- North coordinates include corners
|  |          |  |
|  |          |  | <- East and west coordinates ALSO include corners
|  |          |  |
```

For ease of use, each boundary should have a constructor option where the user can provide
an [`AbstractGrid`](@ref) concrete type and define a boundary where `val` is at the edge of the
grid so that thr boundaries form a border around the grid.

## _API_
...
"""
abstract type AbstractBoundary{
    D<:AbstractDirection,
    FT,
}<:AbstractDomainElement{FT} end

# Concrete type of AbstractBoundary - see documentation below
struct OpenBoundary{D, FT}<:AbstractBoundary{D, FT}
    poly::StaticQuadrilateral{FT}
    val::FT
end

"""
    OpenBoundary{D, FT} <: AbstractBoundary{D, FT}

A sub-type of AbstractBoundary that allows a floe to pass out of the domain edge
without any effects on the floe.
"""
function OpenBoundary(::Type{D}, ::Type{FT} = Float64; grid = nothing, poly = nothing, val = nothing) where {D, FT}
    if !isnothing(grid)
        poly, val = _grid_boundary_info(FT, D, grid)
    elseif isnothing(poly) || isnothing(val)
        throw(ArgumentError("To create an OpenBoundary, either provide a grid or a poly AND a val."))
    end
    return OpenBoundary{D, FT}(poly, val)
end

# Concrete type of AbstractBoundary - see documentation below
struct PeriodicBoundary{D, FT}<:AbstractBoundary{D, FT}
    poly::StaticQuadrilateral{FT}
    val::FT
end

"""
    PeriodicBoundary <: AbstractBoundary

A sub-type of AbstractBoundary that moves a floe from one side of the domain to
the opposite side of the domain when it crosses the boundary, bringing the floe
back into the domain.
"""
function PeriodicBoundary(::Type{D}, ::Type{FT} = Float64; grid = nothing, poly = nothing, val = nothing) where {D, FT}
    if !isnothing(grid)
        poly, val = _grid_boundary_info(FT, D, grid)
    elseif isnothing(poly) || isnothing(val)
        throw(ArgumentError("To create an PeriodicBoundary, either provide a grid or a poly AND a val."))
    end
    return PeriodicBoundary{D, FT}(poly, val)
end

# Concrete type of AbstractBoundary - see documentation below
struct CollisionBoundary{D, FT}<:AbstractBoundary{D, FT}
    poly::StaticQuadrilateral{FT}
    val::FT
end

"""
    CollisionBoundary <: AbstractBoundary

A sub-type of AbstractBoundary that stops a floe from exiting the domain by
having the floe collide with the boundary. The boundary acts as an immovable,
unbreakable ice floe in the collision. 
"""
function CollisionBoundary(::Type{D}, ::Type{FT} = Float64; grid = nothing, poly = nothing, val = nothing) where {D, FT}
    if !isnothing(grid)
        poly, val = _grid_boundary_info(FT, D, grid)
    elseif isnothing(poly) || isnothing(val)
        throw(ArgumentError("To create an CollisionBoundary, either provide a grid or a poly AND a val."))
    end
    return CollisionBoundary{D, FT}(poly, val)
end

# Concrete type of AbstractBoundary - see documentation below
mutable struct MovingBoundary{D, FT}<:AbstractBoundary{D, FT}
    poly::StaticQuadrilateral{FT}
    val::FT
    u::FT
    v::FT
end


"""
    MovingBoundary <: AbstractBoundary

A sub-type of AbstractBoundary that creates a floe along the boundary that moves
towards the center of the domain at the given velocity, compressing the ice
within the domain. This subtype is a mutable struct so that the coordinates and
val can be changed as the boundary moves. The u and v velocities are in [m/s].

Note that with a u-velocity, east and west walls move towards the center of the domain,
providing compressive stress, and with a v-velocity, the bounday creates a shear stress by
incorporating the velocity into friction calculations but doesn't actually move. This means
that the boundaries cannot move at an angle, distorting the shape of the domain regardless
the the combination of u and v velocities. The same, but opposite is true for the north
and south walls.
"""
function MovingBoundary(::Type{D}, ::Type{FT} = Float64; u = 0.0, v = 0.0,
    grid = nothing, poly = nothing, val = nothing,
) where {D, FT}
    if !isnothing(grid)
        poly, val = _grid_boundary_info(FT, D, grid)
    elseif isnothing(poly) || isnothing(val)
        throw(ArgumentError("To create an MovingBoundary, either provide a grid or a poly AND a val."))
    end
    if u == v == 0
        @warn "MovingBoundary velocities are both zero. Boundary will not move. Use keywords u and v to set velocities."
    end
    return MovingBoundary{D, FT}(poly, val, u, v)
end

# Union of all non-peridic boundary types to use as shorthand for dispatch.
const NonPeriodicBoundary = Union{
    OpenBoundary,
    CollisionBoundary,
    MovingBoundary,
}

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

"""
    periodic_compat(b1, b2)

Checks if two boundaries are compatible as a periodic pair. This is true if they
are both periodic, or if neither are periodic. Otherwise, it is false. 
"""
function periodic_compat(::PeriodicBoundary, ::PeriodicBoundary)
    return true
end

function periodic_compat(::PeriodicBoundary, _)
    return false
end

function periodic_compat(_, ::PeriodicBoundary)
    return false
end

function periodic_compat(_, _)
    return true
end

"""
Domain that holds 4 Boundary elements, forming a rectangle bounding the model
during the simulation, and a list of topography elements.

In order to create a Domain, three conditions need to be met. First, if needs to
be periodically compatible. This means that pairs of opposite boundaries both
need to be periodic if one of them is periodic. Next, the value in the north
boundary must be greater than the south boundary and the value in the east
boundary must be greater than the west in order to form a valid rectangle.

Note: The code depends on the boundaries forming a rectangle oriented along the
cartesian grid. Other shapes/orientations are not supported at this time. 
"""
struct Domain{
    FT<:AbstractFloat,
    NB<:AbstractBoundary{North, FT},
    SB<:AbstractBoundary{South, FT},
    EB<:AbstractBoundary{East, FT},
    WB<:AbstractBoundary{West, FT},
    TT<:StructArray{<:TopographyElement{FT}},
}
    north::NB
    south::SB
    east::EB
    west::WB
    topography::TT

    function Domain{FT, NB, SB, EB, WB, TT}(
        north::NB,
        south::SB,
        east::EB,
        west::WB,
        topography::TT,
    ) where {
        FT<:AbstractFloat,
        NB<:AbstractBoundary{North, FT},
        SB<:AbstractBoundary{South, FT},
        EB<:AbstractBoundary{East, FT},
        WB<:AbstractBoundary{West, FT},
        TT<:StructArray{<:TopographyElement{FT}},
    }
        if !periodic_compat(north, south)
            throw(ArgumentError("North and south boundary walls are not \
                periodically compatable as only one of them is periodic."))
        elseif !periodic_compat(east, west)
            throw(ArgumentError("East and west boundary walls are not \
                periodically compatable as only one of them is periodic."))
        elseif north.val < south.val
            throw(ArgumentError("North boundary value is less than south \
                boundary value."))
        elseif east.val < west.val
            throw(ArgumentError("East boundary value is less than west \
                boundary value."))
        end
        new{FT, NB, SB, EB, WB, TT}(north, south, east, west, topography)
    end

    Domain(
        north::NB,
        south::SB,
        east::EB,
        west::WB,
        topography::TT,
    ) where {
        FT<:AbstractFloat,
        NB<:AbstractBoundary{North, FT},
        SB<:AbstractBoundary{South, FT},
        EB<:AbstractBoundary{East, FT},
        WB<:AbstractBoundary{West, FT},
        TT<:StructArray{<:TopographyElement{FT}},
    } =
        Domain{FT, NB, SB, EB, WB, TT}(north, south, east, west, topography)
end

function get_domain_element(domain, idx)
    if idx == -1
        return domain.north
    elseif idx == -2
        return domain.south
    elseif idx == -3
        return domain.east
    elseif idx == -4
        return domain.west
    else
        topo_idx = -(idx + 4)
        return get_floe(domain.topography, topo_idx)
    end
end

"""
    Domain(north, south, east, west)

Creates domain with empty list of topography and given boundaries.
Inputs:
    north   <AbstractBoundary> north boundary
    south   <AbstractBoundary> south boundary
    east    <AbstractBoundary> east boundary
    west    <AbstractBoundary> west boundary
"""
Domain(
    north::NB,
    south::SB,
    east::EB,
    west::WB,
) where {
    FT<:AbstractFloat,
    NB<:AbstractBoundary{North, FT},
    SB<:AbstractBoundary{South, FT},
    EB<:AbstractBoundary{East, FT},
    WB<:AbstractBoundary{West, FT},
} =
    Domain(
        north,
        south,
        east,
        west,
        StructArray{TopographyElement{FT}}(undef, 0),
    )
