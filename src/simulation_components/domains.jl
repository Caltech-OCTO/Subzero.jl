"""
Structs and functions used to define a Subzero domains
"""

"""
    AbstractDirection

An abstract type for the boundary cardinal directions within model domain.
Boundary direction will control behavior of sea ice floes at edges of domain.
"""
abstract type AbstractDirection end

"""
    North<:AbstractDirection

A simple direction type representing if a boundary is the northern boundary in a
rectangular domain.
"""
struct North<:AbstractDirection end

"""
    South<:AbstractDirection

A simple direction type representing if a boundary is the southern boundary in a
rectangular domain.
"""
struct South<:AbstractDirection end

"""
    East<:AbstractDirection

A simple direction type representing if a boundary is the eastern boundary in a
rectangular domain.
"""
struct East<:AbstractDirection end

"""
    West<:AbstractDirection

A simple direction type representing if a boundary is the western boundary in a
rectangular domain.
"""
struct West<:AbstractDirection end

"""
    boundary_coords(grid::AbstractGrid, ::Type{North})

Determine coordinates of northen-most boundary of domain if around the edge of
the grid.
Inputs:
    grid    <AbstractGrid> model grid
            <Type{North}> boundary direction type
Output:
    PolyVec of boundary coordinates. These coordinates describe a rectangle that
    has a length 2-times the length of the grid in the x-direction, centered on
    the grid so that there is a buffer of half of the grid on either side. The
    height is half of the grid in the y-direction. This buffer prevents pieces
    of floes from passing outside the boundary before the next timestep -
    possibly too cautious. If boundary_coords methods are used for each
    direction, corners will be shared between adjacent boundaries. 
"""
function boundary_coords(grid::AbstractGrid, ::Type{North})
    Δx = (grid.xf - grid.x0)/2 # Half of the grid in x
    Δy = (grid.yf - grid.y0)/2 # Half of the grid in y
    return grid.yf,  # val
        [[[grid.x0 - Δx, grid.yf],  # coords
          [grid.x0 - Δx, grid.yf + Δy],
          [grid.xf + Δx, grid.yf + Δy], 
          [grid.xf + Δx, grid.yf], 
          [grid.x0 - Δx, grid.yf]]]
end

"""
    boundary_coords(grid::AbstractGrid, ::Type{South})

Determine coordinates of southern-most boundary of domain if around the edge of
the grid.
Inputs:
    grid    <AbstractGrid> model grid
            <Type{South}> boundary direction type
Output:
    PolyVec of boundary coordinates. See documentation of North method of this
    function for more details. 
"""
function boundary_coords(grid::AbstractGrid, ::Type{South})
    Δx = (grid.xf - grid.x0)/2 # Half of the grid in x
    Δy = (grid.yf - grid.y0)/2 # Half of the grid in y
    return grid.y0,  # val
        [[[grid.x0 - Δx, grid.y0 - Δy],  # coords
          [grid.x0 - Δx, grid.y0],
          [grid.xf + Δx, grid.y0], 
          [grid.xf + Δx, grid.y0 - Δy], 
          [grid.x0 - Δx, grid.y0 - Δy]]]
end

"""
    boundary_coords(grid::AbstractGrid, ::Type{East})

Determine coordinates of eastern-most boundary of domain if around the edge of
the grid.
Inputs:
    grid    <AbstractGrid> model grid
            <Type{East}> boundary direction type
Output:
    PolyVec of boundary coordinates. See documentation of North method of this
    function for more details. 
"""
function boundary_coords(grid::AbstractGrid, ::Type{East})
    Δx = (grid.xf - grid.x0)/2 # Half of the grid in x
    Δy = (grid.yf - grid.y0)/2 # Half of the grid in y
    return grid.xf,  # val
        [[[grid.xf, grid.y0 - Δy],  # coords
          [grid.xf, grid.yf + Δy],
          [grid.xf + Δx, grid.yf + Δy], 
          [grid.xf + Δx, grid.y0 - Δy], 
          [grid.xf, grid.y0 - Δy]]]
end

"""
    boundary_coords(grid::AbstractGrid, ::Type{West})

Determine coordinates of western-most boundary of domain if around the edge of
the grid.
Inputs:
    grid    <AbstractGrid> model grid
            <Type{West}> boundary direction
Output:
    PolyVec of boundary coordinates. See documentation of North method of this
    function for more details. 
"""
function boundary_coords(grid::AbstractGrid, ::Type{West})
    Δx = (grid.xf - grid.x0)/2 # Half of the grid in x
    Δy = (grid.yf - grid.y0)/2 # Half of the grid in y
    return grid.x0,  # val
        [[[grid.x0 - Δx, grid.y0 - Δy],  # coords
          [grid.x0 - Δx, grid.yf + Δy],
          [grid.x0, grid.yf + Δy], 
          [grid.x0, grid.y0 - Δy], 
          [grid.x0 - Δx, grid.y0 - Δy]]]
end

"""
    AbstractDomainElement{FT<:AbstractFloat}

An abstract type for all of the element that create the shape of the domain:
the 4 boundary walls that make up the rectangular domain and the topography
within the domain. 
"""
abstract type AbstractDomainElement{FT<:AbstractFloat} end


"""
    AbstractBoundary{D<:AbstractDirection, FT}<:AbstractDomainElement{FT}

An abstract type for the types of boundaries at the edges of the model domain.
Boundary types will control behavior of sea ice floes at edges of domain.
The direction given by type D denotes which edge of a domain this boundary could
be and type FT is the simulation float type (e.g. Float64 or Float32).

Each boundary type has the coordinates of the boudnary as a field. These should
be shapes that completely seal the domain, and should overlap on the corners as
seen in the example below:
 ________________
|__|____val___|__| <- North coordinates include corners
|  |          |  |
|  |          |  | <- East and west coordinates ALSO include corners
|  |          |  |
Each bounday type also has a field called "val" that holds value that defines
the line y = val or x = val (depending on boundary direction), such that if the
floe crosses that line it would be partially within the boundary. 
"""
abstract type AbstractBoundary{
    D<:AbstractDirection,
    FT,
}<:AbstractDomainElement{FT} end

"""
    OpenBoundary <: AbstractBoundary

A sub-type of AbstractBoundary that allows a floe to pass out of the domain edge
without any effects on the floe.
"""
struct OpenBoundary{D, FT}<:AbstractBoundary{D, FT}
    coords::PolyVec{FT}
    val::FT
end

"""
    OpenBoundary(::Type{FT}, ::Type{D}, args...)

A float type FT can be provided as the first argument of any OpenBoundary
constructor. The second argument D is the directional type of the boundary.
An OpenBoundary of type FT and D will be created by passing all other arguments
to the correct constructor. 
"""
OpenBoundary(::Type{FT}, ::Type{D}, args...) where {
    FT <: AbstractFloat,
    D <: AbstractDirection,
} = OpenBoundary{D, FT}(args...)

"""
    OpenBoundary(::Type{D}, args...)

If a float type isn't specified, OpenBoundary will be of type Float64 and the
correct constructor will be called with all other arguments.
"""
OpenBoundary(::Type{D}, args...) where {D <: AbstractDirection} =
    OpenBoundary{D, Float64}(args...)

"""
    OpenBoundary{D, FT}(grid::AbstractGrid)

Creates open boundary on the edge of the grid, and with the direction as a type.
Edge is determined by direction.
Inputs:
    grid        <AbstractGrid> model grid
Outputs:
    Open Boundary on edge of grid given by direction and type. 
"""
function OpenBoundary{D, FT}(
    grid::AbstractGrid,
) where {D <: AbstractDirection, FT <: AbstractFloat}
    val, coords = boundary_coords(grid, D)
    OpenBoundary{D, FT}(
        coords,
        val,
    )
end

"""
    PeriodicBoundary <: AbstractBoundary

A sub-type of AbstractBoundary that moves a floe from one side of the domain to
the opposite side of the domain when it crosses the boundary, bringing the floe
back into the domain.
"""
struct PeriodicBoundary{D, FT}<:AbstractBoundary{D, FT}
    coords::PolyVec{FT}
    val::FT
end

"""
    PeriodicBoundary(::Type{FT}, ::Type{D}, args...)

A float type FT can be provided as the first argument of any PeriodicBoundary
constructor. The second argument D is the directional type of the boundary.
A PeriodicBoundary of type FT and D will be created by passing all other
arguments to the correct constructor. 
"""
PeriodicBoundary(::Type{FT}, ::Type{D}, args...) where {
    FT <: AbstractFloat,
    D <: AbstractDirection,
} = PeriodicBoundary{D, FT}(args...)

"""
    PeriodicBoundary(::Type{D}, args...)

If a float type isn't specified, PeriodicBoundary will be of type Float64 and
the correct constructor will be called with all other arguments.
"""
PeriodicBoundary(::Type{D}, args...) where {D <: AbstractDirection} =
    PeriodicBoundary{D, Float64}(args...)

"""
    PeriodicBoundary{D, FT}(grid::AbstractGrid)

Creates periodic boundary on the edge of the grid, and with the direction as a
type. Edge is determined by direction.
Inputs:
    grid        <AbstractGrid> model grid
Outputs:
    Periodic Boundary on edge of grid given by direction and type. 
"""
function PeriodicBoundary{D, FT}(
    grid::AbstractGrid,
) where {D <: AbstractDirection, FT <: AbstractFloat}
    val, coords = boundary_coords(grid, D)
    PeriodicBoundary{D, FT}(
        coords,
        val,
    )
end

"""
    CollisionBoundary <: AbstractBoundary

A sub-type of AbstractBoundary that stops a floe from exiting the domain by
having the floe collide with the boundary. The boundary acts as an immovable,
unbreakable ice floe in the collision. 
"""
struct CollisionBoundary{D, FT}<:AbstractBoundary{D, FT}
    coords::PolyVec{FT}
    val::FT
end

"""
    CollisionBoundary(::Type{FT}, ::Type{D}, args...)

A float type FT can be provided as the first argument of any CollisionBoundary
constructor. The second argument D is the directional type of the boundary.
A CollisionBoundary of type FT and D will be created by passing all other
arguments to the correct constructor. 
"""
CollisionBoundary(::Type{FT}, ::Type{D}, args...) where {
    FT <: AbstractFloat,
    D <: AbstractDirection,
} = CollisionBoundary{D, FT}(args...)

"""
    CollisionBoundary(::Type{D}, args...)

If a float type isn't specified, CollisionBoundary will be of type Float64 and
the correct constructor will be called with all other arguments.
"""
CollisionBoundary(::Type{D}, args...) where {D <: AbstractDirection} =
    CollisionBoundary{D, Float64}(args...)


"""
    CollisionBoundary{D, FT}(grid)

Creates collision boundary on the edge of the grid, and with the direction as a
type. Edge is determined by direction.
Inputs:
    grid        <AbstractGrid> model grid
    direction   <AbstractDirection> direction of boundary wall
Outputs:
    Collision Boundary on edge of grid given by direction. 
"""
function CollisionBoundary{D, FT}(
    grid::AbstractGrid,
) where {D <: AbstractDirection, FT <: AbstractFloat}
    val, coords = boundary_coords(grid, D)
    CollisionBoundary{D, FT}(
        coords,
        val,
    )
end

"""
    CompressionBoundary <: AbstractBoundary

A sub-type of AbstractBoundary that creates a floe along the boundary that moves
towards the center of the domain at the given velocity, compressing the ice
within the domain. This subtype is a mutable struct so that the coordinates and
val can be changed as the boundary moves. The velocity is in [m/s].
"""
mutable struct CompressionBoundary{D, FT}<:AbstractBoundary{D, FT}
    coords::PolyVec{FT}
    val::FT
    velocity::FT
end

"""
    CompressionBoundary(::Type{FT}, ::Type{D}, args...)

A float type FT can be provided as the first argument of any CompressionBoundary
constructor. The second argument D is the directional type of the boundary.
A CompressionBoundary of type FT and D will be created by passing all other
arguments to the correct constructor. 
"""
CompressionBoundary(::Type{FT}, ::Type{D}, args...) where {
    FT <: AbstractFloat,
    D <: AbstractDirection,
} = CompressionBoundary{D, FT}(args...)

"""
    CompressionBoundary(::Type{D}, args...)

If a float type isn't specified, CompressionBoundary will be of type Float64 and
the correct constructor will be called with all other arguments.
"""
CompressionBoundary(::Type{D}, args...) where {D <: AbstractDirection} =
    CompressionBoundary{D, Float64}(args...)


"""
    CompressionBoundary{D, FT}(grid, velocity)

Creates compression boundary on the edge of the grid, and with the direction as
a type.
Edge is determined by direction.
Inputs:
        grid        <AbstractGrid> model grid
        velocity    <AbstractFloat> velocity of boundary
Outputs:
    CompressionBoundary on edge of grid given by direction. 
Note:
    If a boundary should move from north to south or east to west, the velocity
    should be negative. Else the velocity should be positive.
"""
function CompressionBoundary{D, FT}(
    grid::AbstractGrid,
    velocity,
) where {D <: AbstractDirection, FT <: AbstractFloat}
    val, coords = boundary_coords(grid, D)
    CompressionBoundary{D, FT}(
        coords,
        val,
        velocity,
    )
end

"""
    NonPeriodicBoundary

Union of all non-peridic boundary types to use as shorthand for dispatch.
"""
const NonPeriodicBoundary = Union{
    OpenBoundary,
    CollisionBoundary,
    CompressionBoundary,
}

"""
    TopographyE{FT}<:AbstractDomainElement{FT}

Singular topographic element with coordinates field storing where the element is
within the grid. These are used to create the desired topography within the
simulation and will be treated as islands or coastline within the model
in that they will not move or break due to floe interactions, but they will
affect floes.
"""
struct TopographyElement{FT}<:AbstractDomainElement{FT}
    coords::PolyVec{FT}
    centroid::Vector{FT}
    rmax::FT

    function TopographyElement{FT}(
        coords,
        centroid,
        rmax,
    ) where {FT <: AbstractFloat}
        if rmax <= 0
            throw(ArgumentError("Topography element maximum radius must be \
                positive and non-zero."))
        end
        new{FT}(valid_polyvec!(rmholes(coords)), centroid, rmax)
    end
end

"""
    TopographyElement(::Type{FT}, args...)

A float type FT can be provided as the first argument of any TopographyElement
constructor. A TopographyElement of type FT will be created by passing all other
arguments to the correct constructor. 
"""
TopographyElement(::Type{FT}, args...) where {FT <: AbstractFloat} =
    TopographyElement{FT}(args...)

"""
    TopographyElement(args...)

If a type isn't specified, TopographyElement will be of type Float64 and the
correct constructor will be called with all other arguments.
"""
TopographyElement(args...) = TopographyElement{Float64}(args...)

"""
    TopographyElement{FT}(poly)

Constructor for topographic element with LibGEOS Polygon
Inputs:
    poly    <LibGEOS.Polygon>
Output:
    Topographic element of abstract float type FT
"""
function TopographyElement{FT}(poly::LG.Polygon) where {FT <: AbstractFloat}
    topo = rmholes(poly)
    centroid = find_poly_centroid(topo)
    coords = find_poly_coords(topo)
    # Move coordinates to be centered at origin to calculate maximum radius
    translate!(
        coords,
        -centroid[1],
        -centroid[2],
    )
    rmax = sqrt(maximum([sum(c.^2) for c in coords[1]]))
    translate!(
        coords,
        centroid[1],
        centroid[2],
    )
    return TopographyElement{FT}(
        coords,
        centroid,
        rmax,
    )
end

"""
    TopographyElement{FT}(coords)

Constructor for topographic element with PolyVec coordinates
Inputs:
    coords      <PolyVec>
Output:
    Topographic element of abstract float type FT
"""
function TopographyElement{FT}(
    coords::PolyVec,
) where {FT <: AbstractFloat} 
    return TopographyElement{FT}(
        LG.Polygon(
            convert(PolyVec{Float64}, coords),
        ),
    )
end

"""
    initialize_topography_field(args...)

If a type isn't specified, the list of TopographyElements will each be of type
Float64 and the correct constructor will be called with all other arguments.
"""
initialize_topography_field(args...) =
    initialize_topography_field(Float64, args...)

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
    ::Type{FT},
    coords,
) where {FT <: AbstractFloat}
    topo_arr = StructArray{TopographyElement{FT}}(undef, length(coords))
    for i in eachindex(coords)
        c = Subzero.rmholes(coords[i])
        topo_arr[i] = TopographyElement{FT}(c)
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
}
    north::NB
    south::SB
    east::EB
    west::WB
    topography::StructArray{TopographyElement{FT}}

    function Domain{FT, NB, SB, EB, WB}(
        north::NB,
        south::SB,
        east::EB,
        west::WB,
        topography::StructArray{TopographyElement{FT}},
    ) where {
        FT<:AbstractFloat,
        NB<:AbstractBoundary{North, FT},
        SB<:AbstractBoundary{South, FT},
        EB<:AbstractBoundary{East, FT},
        WB<:AbstractBoundary{West, FT},
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
        new{FT, NB, SB, EB, WB}(north, south, east, west, topography)
    end

    Domain(
        north::NB,
        south::SB,
        east::EB,
        west::WB,
        topography::StructArray{TopographyElement{FT}},
    ) where {
        FT<:AbstractFloat,
        NB<:AbstractBoundary{North, FT},
        SB<:AbstractBoundary{South, FT},
        EB<:AbstractBoundary{East, FT},
        WB<:AbstractBoundary{West, FT},
    } =
        Domain{FT, NB, SB, EB, WB}(north, south, east, west, topography)
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
    Domain{FT, NB, SB, EB, WB}(
        north,
        south,
        east,
        west,
        StructArray{TopographyElement{FT}}(undef, 0),
    )
