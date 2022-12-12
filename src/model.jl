"""
Structs and functions used to define a Subzero model
"""

"""
    Grid{FT<:AbstractFloat}

Grid splitting the model into distinct rectanglular grid cells where xg are the grid lines in the x-direction (1xn vector)
and yg are the grid lines in the y-direction (mx1 vector). xc and cy are the mid-lines on grid cells in the x and y-direction.
These have dimensions n-1x1 and m-1x1 respectively. The dimension field holds the number of rows (m-1) and columns (n-1) in the grid.
The dimensions of each of the fields must match according to the above definitions.

This struct is also used to create a coarse grid of the model domain.
Ocean and floe data is averaged over this coarse grid and then saved as model output.
"""
struct Grid{FT<:AbstractFloat}
    dims::Tuple{Int, Int}
    xg::Vector{FT}
    yg::Vector{FT}
    xc::Vector{FT}
    yc::Vector{FT}

    Grid(dims, xg, yg, xc, yc) =
        (length(xg) == dims[2]+1 && length(yg) == dims[1]+1 &&
         length(xc) == dims[2] && length(yc) == dims[1]) ?
        new{eltype(xg)}(dims, xg, yg, xc, yc) :
        throw(ArgumentError("Dimension field doesn't match grid dimensions."))
end

"""
    Grid(lx, ux, ly, uy, Δx, Δy, t::Type{T} = Float64)

Construct a rectanglular grid for model given upper and lower bounds for x and y and grid cell dimensions.
Inputs: 
        lx       <Real> lower bound of grid x-direction
        ux       <Real> upper bound of grid x-direction
        ly       <Real> lower bound of grid y-direction
        uy       <Real> upper bound of grid y-direction
        Δx       <Real> length/height of grid cells in x-direction
        Δy       <Real> length/height of grid cells in y-direction
        t        <Type> datatype to convert grid fields - must be a Float!
Output: 
        Grid from lx to ux and height from ly to uy with grid squares of size Δx by Δy
Warning: If Δx doesn't evenly divide x length (lu-lx) or Δy doesn't evenly 
         divide y length (uy-ly) you won't get full size grid. The grid will be "trimmed" to the nearest
         full grid square in both directions.
"""
function Grid(lx, ux, ly, uy, Δx, Δy, t::Type{T} = Float64) where T
    xg = collect(T, lx:Δx:ux) 
    yg = collect(T, ly:Δy:uy)
    nx = length(xg) - 1
    ny = length(yg) - 1
    xc = collect(xg[1]+Δx/2:Δx:xg[end]-Δx/2)
    yc = collect(yg[1]+Δy/2:Δy:yg[end]-Δy/2)
    return Grid((ny, nx), xg, yg, xc, yc)
end


"""
    Grid(Lx, Ly, Δx, Δy, t::Type{T} = Float64)

Construct a rectanglular grid for the model given Lx, Ly and cell dimensions.
Inputs: 
        Lx       <Real> grid length will range from 0 to Lx
        Ly       <Real> grid height will range from y to Ly
        Δx       <Real> length/height of grid cells in x-direction
        Δy       <Real> length/height of grid cells in y-direction
        t        <Type> datatype to convert grid fields - must be a Float!
Output: 
        Grid with length of Lx (0.0 to LX) and height of Ly (0.0 to LY) with square Δx by Δy grid cells
Warning: If Δx doesn't evenly divide Lx or Δy doesn't evenly divide Ly you 
         won't get full size grid. The grid will be "trimmed" to the nearest full grid square.
"""
function Grid(Lx, Ly, Δx, Δy, t::Type{T} = Float64) where T
    return Grid(0.0, Lx, 0.0, Ly, Δx, Δy, T)
end

"""
    Grid(Lx, Ly, Δx, Δy, t::Type{T} = Float64)

Construct a rectanglular grid for model given upper and lower bounds for x and y and the number of grid cells
in both the x and y direction.
Inputs: 
        lx       <Real> lower bound of grid x-direction
        ux       <Real> upper bound of grid x-direction
        ly       <Real> lower bound of grid y-direction
        uy       <Real> upper bound of grid y-direction
        dims     <(Int, Int)> grid dimensions - rows -> ny, cols -> nx
        t        <Type> datatype to convert grid fields - must be a Float!
Output: 
        Grid from lx to ux and height from ly to uy with nx grid cells in the x-direction and ny grid cells in the y-direction.
"""
function Grid(lx, ux, ly, uy, dims::Tuple{Int, Int}, t::Type{T} = Float64) where T
    Δx = (ux-lx)/dims[2]
    Δy = (uy-ly)/dims[1]
    return Grid(lx, ux, ly, uy, Δx, Δy, T)
end

"""
    Grid(Lx, Ly, nx, ny, t::Type{T} = Float64)

    Construct a rectanglular grid for the model given Lx, Ly, and the number of grid cells in both the x and y direction.
Inputs: 
        Lx       <Real> grid length will range from 0 to Lx
        Ly       <Real> grid height will range from y to Ly
        dims     <(Int, Int)> grid dimensions
        t        <Type> datatype to convert grid fields - must be a Float!
Output: 
        Grid from 0 to Lx and height from 0 to Ly with nx grid cells in the x-direction and ny grid cells in the y-direction.
"""
function Grid(Lx, Ly, dims::Tuple{Int, Int}, t::Type{T} = Float64) where T
    Δx = Lx/dims[2]
    Δy = Ly/dims[1]
    return Grid(0.0, Lx, 0.0, Ly, Δx, Δy, T)
end

"""
Ocean velocities in the x-direction (u) and y-direction (v). u and v should match the size of the corresponding
model grid so that there is one x and y velocity value for each grid cell. Ocean also needs temperature at the
ocean/ice interface in each grid cell. Ocean fields must all be matricies with the same dimensions.
Model cannot be constructed if size of ocean and grid do not match.
"""
struct Ocean{FT<:AbstractFloat}
    u::Matrix{FT}
    v::Matrix{FT}
    temp::Matrix{FT}
    hflx::Matrix{FT} 
    τx::Matrix{FT}
    τy::Matrix{FT}
    si_area::Matrix{FT}

    Ocean(u, v, temp, hflx, τx, τy, si_frac) =
        (size(u) == size(v) == size(temp) == size(hflx) == size(τx) == size(τy) ==
         size(si_frac)) ?
        new{eltype(u)}(u, v, temp, hflx, τx, τy, si_frac) :
        throw(ArgumentError("All ocean fields matricies must have the same dimensions."))
end

"""
    Ocean(grid, u, v, temp, FT)

Construct model ocean.
Inputs: 
        grid    <Grid> model grid cell
        u       <Real> ocean x-velocity for each grid cell
        v       <Real> ocean y-velocity for each grid cell
        temp    <Real> temperature at ocean/ice interface per grid cell
        t       <Type> datatype to convert ocean fields - must be a Float!
Output: 
        Ocean with constant velocity and temperature in each grid cell.
"""
function Ocean(grid::Grid, u, v, temp, ::Type{T} = Float64) where T
    nvals = grid.dims .+ 1  # one value per grid line - not grid cell 
    return Ocean(fill(convert(T, u), nvals), 
                 fill(convert(T, v), nvals), 
                 fill(convert(T, temp), nvals),
                 zeros(T, nvals), zeros(T, nvals), 
                 zeros(T, nvals), zeros(T, nvals))
end

"""
Wind velocities in the x-direction (u) and y-direction (v). u and v should match the size of the corresponding
model grid so that there is one x and y velocity value for each grid cell. Wind also needs temperature at the
atmosphere/ice interface in each grid cell. Model cannot be constructed if size of wind and grid do not match.
"""
struct Wind{FT<:AbstractFloat}
    u::Matrix{FT}
    v::Matrix{FT}
    temp::Matrix{FT}

    Wind(u, v, temp) =
    (size(u) == size(v) == size(temp)) ?
    new{eltype(u)}(u, v, temp) :
    throw(ArgumentError("All wind fields matricies must have the same dimensions."))
end

"""
    Wind(grid, u, v, FT)

Construct model atmosphere/wind.
Inputs: 
        grid    <Grid> model grid cell
        u       <Real> wind x-velocity for each grid cell
        v       <Real> wind y-velocity for each grid cell
        temp    <Real> temperature at atmopshere/ice interface per grid cell
        t       <Type> datatype to convert ocean fields - must be a Float!
Output: 
        Ocean with constant velocity and temperature in each grid cell.
"""
function Wind(grid, u, v, temp, ::Type{T} = Float64) where T
    nvals = grid.dims .+ 1  # one value per grid line - not grid cell 
    return Wind(fill(convert(T, u), nvals),
                fill(convert(T, v), nvals),
                fill(convert(T, temp), nvals))
end

"""
    AbstractDirection

An abstract type for the boundary cardinal directions within model domain.
Boundary direction will control behavior of sea ice floes at edges of domain.
"""
abstract type AbstractDirection end

"""
    North<:AbstractDirection

A simple direction type representing if a boundary is the northern boundary in a rectangular domain.
"""
struct North<:AbstractDirection end

"""
    South<:AbstractDirection

A simple direction type representing if a boundary is the souther boundary in a rectangular domain.
"""
struct South<:AbstractDirection end

"""
    East<:AbstractDirection

A simple direction type representing if a boundary is the eastern boundary in a rectangular domain.
"""
struct East<:AbstractDirection end

"""
    West<:AbstractDirection

A simple direction type representing if a boundary is the western boundary in a rectangular domain.
"""
struct West<:AbstractDirection end

"""
    boundary_coords(grid::Grid, ::North)

Determine coordinates of northen-most boundary of domain if around the edge of the grid.
Inputs:
        grid    <Grid> model grid
                <North> boundary direction
Output:
        PolyVec of boundary coordinates. These coordinates describe a rectangle that has a length 2-times
        the length of the grid in the x-direction, centered on the grid so that there is a buffer of half 
        of the grid on either side. The height is half of the grid in the y-direction. This buffer prevents
        pieces of floes from passing outside the boundary before the next timestep - possibly too cautious.
        If boundary_coords methods are used for each direction, corners will be shared between adjacent boundaries. 
"""
function boundary_coords(grid::Grid, ::North)
    Δx = (grid.xg[end] - grid.xg[1])/2 # Half of the grid in x
    Δy = (grid.yg[end] - grid.yg[1])/2 # Half of the grid in y
    return grid.yg[end],  # val
        [[[grid.xg[1] - Δx, grid.yg[end]],  # coords
          [grid.xg[1] - Δx, grid.yg[end] + Δy],
          [grid.xg[end] + Δx, grid.yg[end] + Δy], 
          [grid.xg[end] + Δx, grid.yg[end]], 
          [grid.xg[1] - Δx, grid.yg[end]]]]
end

"""
    boundary_coords(grid::Grid, ::South)

Determine coordinates of southern-most boundary of domain if around the edge of the grid.
Inputs:
        grid    <Grid> model grid
                <South> boundary direction
Output:
        PolyVec of boundary coordinates. See documentation of North method of this function for more details. 
"""
function boundary_coords(grid::Grid, ::South)
    Δx = (grid.xg[end] - grid.xg[1])/2 # Half of the grid in x
    Δy = (grid.yg[end] - grid.yg[1])/2 # Half of the grid in y
    return grid.yg[1],  # val
        [[[grid.xg[1] - Δx, grid.yg[1] - Δy],  # coords
          [grid.xg[1] - Δx, grid.yg[1]],
          [grid.xg[end] + Δx, grid.yg[1]], 
          [grid.xg[end] + Δx, grid.yg[1] - Δy], 
          [grid.xg[1] - Δx, grid.yg[1] - Δy]]]
end

"""
    boundary_coords(grid::Grid, ::East)

Determine coordinates of eastern-most boundary of domain if around the edge of the grid.
Inputs:
        grid    <Grid> model grid
                <East> boundary direction
Output:
        PolyVec of boundary coordinates. See documentation of North method of this function for more details. 
"""
function boundary_coords(grid::Grid, ::East)
    Δx = (grid.xg[end] - grid.xg[1])/2 # Half of the grid in x
    Δy = (grid.yg[end] - grid.yg[1])/2 # Half of the grid in y
    return grid.xg[end],  # val
        [[[grid.xg[end], grid.yg[1] - Δy],  # coords
          [grid.xg[end], grid.yg[end] + Δy],
          [grid.xg[end] + Δx, grid.yg[end] + Δy], 
          [grid.xg[end] + Δx, grid.yg[1] - Δy], 
          [grid.xg[end], grid.yg[1] - Δy]]]
end

"""
    boundary_coords(grid::Grid, ::West)

Determine coordinates of western-most boundary of domain if around the edge of the grid.
Inputs:
        grid    <Grid> model grid
                <West> boundary direction
Output:
        PolyVec of boundary coordinates. See documentation of North method of this function for more details. 
"""
function boundary_coords(grid::Grid, ::West)
    Δx = (grid.xg[end] - grid.xg[1])/2 # Half of the grid in x
    Δy = (grid.yg[end] - grid.yg[1])/2 # Half of the grid in y
    return grid.xg[1],  # val
        [[[grid.xg[1] - Δx, grid.yg[1] - Δy],  # coords
          [grid.xg[1] - Δx, grid.yg[end] + Δy],
          [grid.xg[1], grid.yg[end] + Δy], 
          [grid.xg[1], grid.yg[1] - Δy], 
          [grid.xg[1] - Δx, grid.yg[1] - Δy]]]
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
The direction given by type D denotes which edge of a domain this boundary could be.

Each boundary type has the coordinates of the boudnary as a field. These should
be shapes that completely seal the domain, and should overlap on the corners as
seen  in the example below:
 ________________
|__|____val___|__| <- North coordinates must include corners
|  |          |  |
|  |          |  | <- East coordinates must ALSO include corners
|  |          |  |
Each bounday type also has a field called "val" that holds value that defines
the line y = val or x = val (depending on boundary direction), such that if the
floe crosses that line it would be partially within the boundary. 
"""
abstract type AbstractBoundary{D<:AbstractDirection, FT}<:AbstractDomainElement{FT} end

"""
    OpenBoundary <: AbstractBoundary

A sub-type of AbstractBoundary that allows a floe to pass out of the domain edge without any effects on the floe.
"""
struct OpenBoundary{D, FT}<:AbstractBoundary{D, FT}
    coords::PolyVec{FT}
    val::FT
end

"""
    OpenBoundary(coords, val, direction::AbstractDirection)

Creates open boundary with given coords and val, and with the direction as a type.
Inputs:
        coords      <PolyVec{AbstractFloat}> coordinates of boundary
        val         <AbstractFloat> value defining line that marks edge of domain
        direction   <AbstractDirection> direction of boundary wall
Output:
        Open Boundary of type given by direction and defined by given coordinates and edge value
"""
function OpenBoundary(coords, val, direction::AbstractDirection)
    return OpenBoundary{typeof(direction), typeof(val)}(coords, val)
end

"""
    OpenBoundary(grid::Grid, direction)

Creates open boundary on the edge of the grid, and with the direction as a type.
Edge is determined by direction.
Inputs:
        grid        <Grid> model grid
        direction   <AbstractDirection> direction of boundary wall
Outputs:
        Open Boundary on edge of grid given by direction. 
"""
function OpenBoundary(grid::Grid, direction)
    val, coords = boundary_coords(grid, direction)
    OpenBoundary(coords, val, direction)
end

"""
    PeriodicBoundary <: AbstractBoundary

A sub-type of AbstractBoundary that moves a floe from one side of the domain to the opposite
side of the domain when it crosses the boundary, bringing the floe back into the domain.

NOTE: Not implemented yet!
"""
struct PeriodicBoundary{D, FT}<:AbstractBoundary{D, FT}
    coords::PolyVec{FT}
    val::FT
end

"""
    PeriodicBoundary(coords, val, direction::AbstractDirection)

Creates periodic boundary with given coords and val, and with the direction as a type.
Inputs:
        coords      <PolyVec{AbstractFloat}> coordinates of boundary
        val         <AbstractFloat> value defining line that marks edge of domain
        direction   <AbstractDirection> direction of boundary wall
Output:
        Periodic Boundary of type given by direction and defined by given coordinates and edge value
"""
function PeriodicBoundary(coords, val, direction::AbstractDirection)
    return PeriodicBoundary{typeof(direction), typeof(val)}(coords, val)
end

"""
    PeriodicBoundary(grid::Grid, direction)

Creates periodic boundary on the edge of the grid, and with the direction as a type.
Edge is determined by direction.
Inputs:
        grid        <Grid> model grid
        direction   <AbstractDirection> direction of boundary wall
Outputs:
        Periodic Boundary on edge of grid given by direction. 
"""
function PeriodicBoundary(grid::Grid, direction)
    val, coords = boundary_coords(grid, direction)
    PeriodicBoundary(coords, val, direction)
end

"""
    CollisionBoundary <: AbstractBoundary

A sub-type of AbstractBoundary that stops a floe from exiting the domain by having the floe collide with the boundary.
The boundary acts as an immovable, unbreakable ice floe in the collision. 
"""
struct CollisionBoundary{D, FT}<:AbstractBoundary{D, FT}
    coords::PolyVec{FT}
    val::FT
end

"""
    CollisionBoundary(coords, val, direction::AbstractDirection)

Creates collision boundary with given coords and val, and with the direction as a type.
Inputs:
        coords      <PolyVec{AbstractFloat}> coordinates of boundary
        val         <AbstractFloat> value defining line that marks edge of domain
        direction   <AbstractDirection> direction of boundary wall
Output:
        Collision Boundary of type given by direction and defined by given coordinates and edge value
"""
function CollisionBoundary(coords, val, direction::AbstractDirection)
    return CollisionBoundary{typeof(direction), typeof(val)}(coords, val)
end

"""
    CollisionBoundary(grid::Grid, direction)

Creates collision boundary on the edge of the grid, and with the direction as a type.
Edge is determined by direction.
Inputs:
        grid        <Grid> model grid
        direction   <AbstractDirection> direction of boundary wall
Outputs:
        Collision Boundary on edge of grid given by direction. 
"""
function CollisionBoundary(grid::Grid, direction)
    val, coords = boundary_coords(grid, direction)
    CollisionBoundary(coords, val, direction)
end
"""
    CompressionBC <: AbstractBC

A sub-type of AbstractBoundary that creates a floe along the boundary that moves towards
the center of the domain at the given velocity, compressing the ice within the domain. This subtype is
a mutable struct so that the coordinates and val can be changed as the boundary moves.
The velocity is in [m/s].

NOTE: Not implemented yet!
"""
mutable struct CompressionBoundary{D, FT}<:AbstractBoundary{D, FT}
    coords::PolyVec{FT}
    val::FT
    velocity::FT
end

"""
    CompressionBoundary(coords, val, direction::AbstractDirection)

Creates compression boundary with given coords and val, and with the direction as a type.
Inputs:
        coords      <PolyVec{AbstractFloat}> coordinates of boundary
        val         <AbstractFloat> value defining line that marks edge of domain
        direction   <AbstractDirection> direction of boundary wall
Output:
        Compression Boundary of type given by direction and defined by given coordinates and edge value
"""
function CompressionBoundary(coords, val, velocity, direction::AbstractDirection)
    return CompressionBoundary{typeof(direction), typeof(val)}(coords, val, velocity)
end

"""
    CompressionBoundary(grid::Grid, direction)

Creates compression boundary on the edge of the grid, and with the direction as a type.
Edge is determined by direction.
Inputs:
        grid        <Grid> model grid
        direction   <AbstractDirection> direction of boundary wall
Outputs:
        Open Boundary on edge of grid given by direction. 
"""
function CompressionBoundary(grid::Grid, direction, velocity)
    val, coords = boundary_coords(grid, direction)
    CompressionBoundary(coords, val, velocity, direction)
end


"""
    TopographyE{FT}<:AbstractDomainElement{FT}

Singular topographic element with coordinates field storing where the element is
within the grid. These are used to create the desired topography within the
simulation and will be treated as islands or coastline within the model
in that they will not move or break due to floe interactions, but they will affect floes.
"""
struct TopographyElement{FT}<:AbstractDomainElement{FT}
    coords::PolyVec{FT}
    centroid::Vector{FT}
    rmax::FT

    function TopographyElement{FT}(coords::PolyVec{FT}, centroid::Vector{FT}, rmax::FT) where {FT <: AbstractFloat}
        rmax > 0 ? 
            new{FT}(valid_polyvec!(rmholes(coords)), centroid, rmax) : 
            throw(ArgumentError("Topography element maximum radius must be positive and non-zero."))
    end

    TopographyElement(coords::PolyVec{FT}, centroid::Vector{FT}, rmax::FT) where {FT <: AbstractFloat} =
        TopographyElement{FT}(coords, centroid, rmax)
end

"""
    Topography(poly::LG.Polygon, t::Type{T} = Float64)

Constructor for topographic element with LibGEOS Polygon
    Inputs:
            poly    <LibGEOS.Polygon> 
            _       <Type> datatype used to run model (Float32 or Float64)
    Output:
            Topographic element coordinates of type PolyVec{T} with any holes removed
    Note:
            Types are specified at Float64 below as type annotations given that when written LibGEOS could exclusivley use Float64 (as of 09/29/22). When this is fixed, this annotation will need to be updated.
            We should only run the model with Float64 right now or else we will be converting the Polygon back and forth all of the time. 
"""
function TopographyElement(poly::LG.Polygon, ::Type{T} = Float64) where T
    topo = rmholes(poly)
    centroid = convert(Vector{T}, LG.GeoInterface.coordinates(LG.centroid(topo))::Vector{Float64})
    coords = convert(PolyVec{T}, LG.GeoInterface.coordinates(topo)::PolyVec{Float64})
    origin_coords = translate(coords, -centroid)
    rmax = sqrt(maximum([sum(c.^2) for c in origin_coords[1]]))
    return TopographyElement(coords, centroid, rmax)
end

"""
    Topography(coords::PolyVec{T}, ::Type{T} = Float64)

Constructor for topographic element with PolyVec coordinates
    Inputs:
            poly    <LibGEOS.Polygon> 
            _       <Type> datatype used to run model (Float32 or Float64)
    Output:
            Topographic element coordinates of type PolyVec{T} with any holes removed
"""
function TopographyElement(coords::PolyVec{T}, ::Type{T} = Float64) where {T} 
    return TopographyElement(LG.Polygon(coords), T)
end

"""
    periodic_compat(b1, b2)

Checks if two boundaries are compatible as a periodic pair. This is true if they are both periodic,
or if neither are periodic. Otherwise, it is false. 
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
Domain that holds 4 Boundary elements, forming a rectangle bounding the model during the simulation. 

In order to create a Domain, three conditions need to be met. First, if needs to be periodically compatible.
This means that pairs of opposite boundaries both need to be periodic if one of them is periodic.
Next, the value in the north boundary must be greater than the south boundary and the value in the east boundary
must be greater than the west in order to form a valid rectangle.

Note: The code depends on the boundaries forming a rectangle oriented along the
cartesian grid. Other shapes/orientations are not supported at this time. 
"""
struct Domain{FT<:AbstractFloat, NB<:AbstractBoundary{North, FT}, SB<:AbstractBoundary{South, FT},
EB<:AbstractBoundary{East, FT}, WB<:AbstractBoundary{West, FT}}
    north::NB
    south::SB
    east::EB
    west::WB
    topography::StructArray{TopographyElement{FT}}

    Domain(north, south, east, west, topography) = 
        (periodic_compat(north, south) && periodic_compat(east, west)) &&
        (north.val > south.val && east.val > west.val) ?
        new{typeof(north.val), typeof(north), typeof(south), typeof(east), typeof(west)}(north, south, east, west, topography) : 
        throw(ArgumentError("Periodic boundary must have matching opposite boundary and/or North value must be greater then South and East must be greater than West."))
end


Domain(north, south, east, west, ::Type{T} = Float64) where T =
    Domain(north, south, east, west, StructArray{TopographyElement{T}}(undef, 0, 0))


"""
Singular sea ice floe with fields describing current state.
"""
@kwdef mutable struct Floe{FT<:AbstractFloat}
    # Physical Properties
    centroid::Vector{FT}    # center of mass of floe (might not be in floe!)
    coords::PolyVec{FT}     # floe coordinates
    height::FT              # floe height (m)
    area::FT                # floe area (m^2)
    mass::FT                # floe mass (kg)
    rmax::FT                # distance of vertix farthest from centroid (m)
    moment::FT              # mass moment of intertia
    angles::Vector{FT}      # interior angles of floe
    # Monte Carlo Points
    mc_x::Vector{FT}        # x-coordinates for monte carlo integration
    mc_y::Vector{FT}        # y-coordinates for monte carlo integration
    # Velocity/Orientation
    α::FT = 0.0             # floe rotation from starting position in radians
    u::FT = 0.0             # floe x-velocity
    v::FT = 0.0             # floe y-velocity
    ξ::FT = 0.0             # floe angular velocity
    # Status
    alive::Int = 1          # floe is still active in simulation
    # Forces/Collisions
    fxOA::FT = 0.0          # force from ocean and wind in x direction
    fyOA::FT = 0.0          # force from ocean and wind in y direction
    trqOA::FT = 0.0      # torque from ocean and wind
    hflx::FT = 0.0          # heat flux under the floe
    overarea::FT = 0.0      # total overlap with other floe
    collision_force::Matrix{FT} = [0.0 0.0] 
    collision_trq::FT = 0.0
    interactions::NamedMatrix{FT} = NamedArray(zeros(7),
    (["floeidx", "xforce", "yforce", "xpoint", "ypoint", "torque", "overlap"]))'
    # Previous values for timestepping
    p_dxdt::FT = 0.0        # previous timestep x-velocity
    p_dydt::FT = 0.0        # previous timestep y-velocity
    p_dudt::FT = 0.0        # previous timestep x-acceleration
    p_dvdt::FT = 0.0        # previous timestep x-acceleration
    p_dξdt::FT = 0.0        # previous timestep time derivative of ξ
    p_dαdt::FT = 0.0        # previous timestep ξ
end

function generate_mc_points(npoints, xfloe, yfloe, rmax, area, t::Type{T} = Float64) where T
    mc_x = rmax * (2rand(T, Int(npoints)) .- 1)
    mc_y = rmax * (2rand(T, Int(npoints)) .- 1)
    mc_in = inpoly2(hcat(mc_x, mc_y), hcat(xfloe, yfloe))
    err = abs(sum(mc_in)/npoints * 4 * rmax^2 - area)/area
    return mc_x, mc_y, mc_in, err
end

"""
    Floe(poly::LG.Polygon, hmean, Δh, ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, t::Type{T} = Float64)

Constructor for floe with LibGEOS Polygon
Inputs:
        poly        <LibGEOS.Polygon> 
        h_mean      <Real> mean height for floes
        Δh          <Real> variability in height for floes
        grid        <Grid>
        ρi          <Real> ice density kg/m3 - default 920
        u           <Real> x-velocity of the floe - default 0.0
        v           <Real> y-velcoity of the floe - default 0.0
        ksi         <Real> angular velocity of the floe - default 0.0
        mc_n        <Real> number of monte carlo points
        t           <Type> datatype to convert ocean fields - must be a Float!
Output:
        Floe with needed fields defined - all default field values used so all forcings start at 0 and floe is "alive".
        Velocities and the density of ice can be optionally set.
Note:
        Types are specified at Float64 below as type annotations given that when written LibGEOS could exclusivley use Float64 (as of 09/29/22).
        When this is fixed, this annotation will need to be updated.
        We should only run the model with Float64 right now or else we will be converting the Polygon back and forth all of the time. 
"""
function Floe(poly::LG.Polygon, hmean, Δh; ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, mc_n = 1000.0, t::Type{T} = Float64) where T
    floe = rmholes(poly)
    centroid = LG.GeoInterface.coordinates(LG.centroid(floe))::Vector{Float64}
    h = hmean + (-1)^rand(0:1) * rand() * Δh  # floe height
    area = LG.area(floe)::Float64  # floe area
    mass = area * h * ρi  # floe mass
    moment = calc_moment_inertia(floe, h, ρi = ρi)
    coords = LG.GeoInterface.coordinates(floe)::PolyVec{Float64}
    origin_coords = translate(coords, -centroid)
    ox, oy = seperate_xy(origin_coords)
    rmax = sqrt(maximum([sum(c.^2) for c in origin_coords[1]]))
    angles = calc_poly_angles(coords, T)
    # Generate Monte Carlo Points
    count = 1
    alive = 1
    mc_x, mc_y, mc_in, err = generate_mc_points(mc_n, ox, oy, rmax, area, T)
    while err > 0.1
        mc_x, mc_y, mc_in, err = generate_mc_points(mc_n, ox, oy, rmax, area, T)
        count += 1
        if count > 10
            err = 0.0
            alive = 0
        end
    end
    mc_x = mc_x[mc_in[:, 1] .|  mc_in[:, 2]]
    mc_y = mc_y[mc_in[:, 1] .|  mc_in[:, 2]]

    return Floe(centroid = convert(Vector{T}, centroid), coords = convert(PolyVec{T}, coords),
                height = convert(T, h), area = convert(T, area), mass = convert(T, mass),
                rmax = convert(T, rmax), moment = convert(T, moment), angles = angles,
                u = convert(T, u), v = convert(T, v), ξ = convert(T, ξ),
                mc_x = mc_x, mc_y = mc_y, alive = alive)
end

"""
    Floe(coords::PolyVec{Float64}, h_mean, Δh, ρi = 920.0, u = 0.0,
    v = 0.0, ξ = 0.0, t::Type{T} = Float64) where T

Floe constructor with PolyVec{Float64}(i.e. Vector{Vector{Vector{Float64}}}) coordinates
Inputs:
        coords      <Vector{Vector{Vector{Float64}}}> floe coordinates
        h_mean      <Real> mean height for floes
        Δh          <Real> variability in height for floes
        grid        <Grid>
        rho_ice     <Real> ice density kg/m3
        u           <Real> x-velocity of the floe - default 0.0
        v           <Real> y-velcoity of the floe - default 0.0
        ksi         <Real> angular velocity of the floe - default 0.0
        t           <Type> datatype to convert ocean fields - must be a Float!
Output:
        Floe with needed fields defined - all default field values used so all forcings and velocities start at 0 and floe is "alive"
"""
Floe(coords::PolyVec{<:Real}, h_mean, Δh; ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, mc_n = 1000.0, t::Type{T} = Float64) where T =
    Floe(LG.Polygon(convert(PolyVec{Float64}, coords)), h_mean, Δh, ρi, u, v, ξ, mc_n, T) 
    # Polygon convert is needed since LibGEOS only takes Float64 - when this is fixed convert can be removed

"""
    domain_in_grid(domain::Domain, grid)

Checks if given rectangular domain is within given grid and gives user a warning if domain is not of maximum possible size given grid dimensions.

Inputs:
        domain      <RectangularDomain>
        grid        <Grid>
Outputs:
        <Boolean> true if domain is within grid bounds, else false
"""
function domain_in_grid(domain::Domain, grid)
    northval = domain.north.val
    southval = domain.south.val
    eastval = domain.east.val
    westval = domain.west.val
    if (northval <= grid.yg[end] &&
        southval >= grid.yg[1] &&
        eastval <= grid.xg[end] &&
        westval >= grid.xg[1])
        if (northval != grid.yg[end] ||
            southval != grid.yg[1] ||
            eastval != grid.xg[end] ||
            westval != grid.xg[1])
            @warn "At least one wall of domain is smaller than grid. This could lead to unneeded computation. Consider making grid smaller or domain larger."
        end 
        return true
    end
    return false
end

"""
Model which holds grid, ocean, wind structs, each with the same underlying float type (either Float32 of Float64) and size.
It also holds the domain information, which includes the topography and the boundaries.
Finally, it holds an StructArray of floe structs, again each relying on the same underlying float type.
"""
struct Model{FT<:AbstractFloat, DT<:Domain{FT, <:AbstractBoundary, <:AbstractBoundary, <:AbstractBoundary, <:AbstractBoundary}}
    grid::Grid{FT}
    ocean::Ocean{FT}
    wind::Wind{FT}
    domain::DT
    floes::StructArray{Floe{FT}}

    Model(grid, ocean, wind, domain, floes) =
        ((grid.dims .+ 1) == size(ocean.u) == size(wind.u) && domain_in_grid(domain, grid)) ?
        new{eltype(ocean.u), typeof(domain)}(grid, ocean, wind, domain, floes) :
        throw(ArgumentError("Size of grid does not match size of ocean and/or wind OR domain is not within grid."))
end