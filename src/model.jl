"""
Structs and functions used to define a Subzero model
"""

"""
    Grid{FT<:AbstractFloat}

Grid splitting the model into distinct rectanglular grid cells where xg are the grid lines in the x-direction (1xn vector) and yg are the grid lines in the y-direction (mx1 vector). xc and cy are the mid-lines on grid cells in the x and y-direction. These have dimensions n-1x1 and m-1x1 respectively. The dimension field holds the number of rows (m-1) and columns (n-1) in the grid. The dimensions of each of the fields must match according to the above definitions.

This struct is also used to create a coarse grid of the model domain. Ocean and floe data is averaged over this coarse grid and then saved as model output.
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
         divide y length (uy-ly) you won't get full size grid. The grid will be "trimmed" to the nearest full grid square in both directions.
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

Construct a rectanglular grid for model given upper and lower bounds for x and y and the number of grid cells in both the x and y direction.
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
Ocean velocities in the x-direction (u) and y-direction (v). u and v should match the size of the corresponding model grid so that there is one x and y velocity value for each grid cell. Ocean also needs temperature at the ocean/ice interface in each grid cell. Ocean fields must all be matricies with the same dimensions. Model cannot be constructed if size of ocean and grid do not match.
"""
struct Ocean{FT<:AbstractFloat}
    u::Matrix{FT}
    v::Matrix{FT}
    temp::Matrix{FT}
    τx::Matrix{FT}
    τy::Matrix{FT}
    si_frac::Matrix{FT}

    Ocean(u, v, temp, τx, τy, si_frac) =
        (size(u) == size(v) == size(temp) == size(τx) == size(τy) ==
         size(si_frac)) ?
        new{eltype(u)}(u, v, temp, τx, τy, si_frac) :
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
Ocean(grid::Grid, u, v, temp, t::Type{T} = Float64) where T =
    Ocean(fill(convert(T, u), grid.dims), 
          fill(convert(T, v), grid.dims), 
          fill(convert(T, temp), grid.dims),
          zeros(T, grid.dims), zeros(T, grid.dims), zeros(T, grid.dims))

# TODO: Do we want to be able to use a psi function? - Ask Mukund

"""
Wind velocities in the x-direction (u) and y-direction (v). u and v should match the size of the corresponding model grid so that there is one x and y velocity value for each grid cell. Wind also needs temperature at the atmosphere/ice interface in each grid cell. Model cannot be constructed if size of wind and grid do not match.
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
Wind(grid, u, v, temp, t::Type{T} = Float64) where T = 
    Wind(fill(convert(T, u), grid.dims),
         fill(convert(T, v), grid.dims),
         fill(convert(T, temp), grid.dims))

"""
    AbstractBC

An abstract type for types of boundary conditions for edges of model domain. Boundary conditions will control behavior of sea ice floes at edges of domain.
"""
abstract type AbstractBC end

"""
    OpenBC <: AbstractBC

A simple concrete type of boundary condition, which allows a floe to pass out of the domain edge without any effects on the floe.
"""
struct OpenBC<:AbstractBC end

"""
    PeriodicBC <: AbstractBC

A simple concrete type of boundary condition, which moves a floe from one side of the domain to the opposite side of the domain, bringing the floe back into the grid.

NOTE: Not implemented yet!
"""
struct PeriodicBC<:AbstractBC end 

"""
    CollisionBC <: AbstractBC

A simple concrete type of boundary condition, which stops a floe from exiting the domain. If the COLLISION flag is on within the simulation, this will exert a potentially breaking force on the floe, else it will just stop the floe from leaving the domain.

NOTE: Not implemented yet!
"""
struct CollisionBC<:AbstractBC end

"""
    CompressionBC <: AbstractBC

A simple concrete type of boundary condition, which creates a floe along the boundary that moves from the boundary towards the center of the domain, compressing the ice within the dominan.

NOTE: Not implemented yet!
"""
struct CompressionBC<:AbstractBC end

"""
    Boundary{BC<:AbstractBC, FT<:AbstractFloat}

One straight edge of a rectangular domian. Holds the boundary condition for the edge and a value that represents the constant x or y-value that the boundary is located at. This value is assigned as an x or a y value based on if the boundary is a north/south "wall" of a rectangular domain (y-value) or a east/west "wall" of the domain (x-value). 

For example, if a boundary has a condition of "open" and a value of 5.0 and is then the "north" edge of a rectangular domain it will be the top of the rectangle at y = 5.0 and floes will pass through the top of the domain without being acted on by any forces.
"""
struct Boundary{BC<:AbstractBC, FT<:AbstractFloat}
    bc::BC
    val::FT 
end

"""
    Boundary(bc, val, t::Type{T} = Float64)

Creates a boundary with given boundary condition and value and converts the value to desired Float type.
Inputs:
        bc      <AbstractBC subtype>
        val     <Real>
        t       <Type> datatype to convert ocean fields - must be a Float! 
Output:
        Boundary of type bc with value of given Float type.
"""
function Boundary(bc, val, t::Type{T} = Float64) where T
    return Boundary(bc, convert(T, val))
end
"""
    AbstractDomain{FT <: AbstractFloat}

An abstract type for shapes of domains that contain the model.
"""
abstract type AbstractDomain{FT <: AbstractFloat} end

"""
    periodic_compat(b1::B, b2::B)

Checks if two boundaries with the same boundary condition B are compatible as opposites. They are by definition, given that if both floes are periodic they work as a north/south pair or a east/west pair and otherwise there are no restrictions on "matching" boundaries.
"""
function periodic_compat(b1::B, b2::B) where B<:Boundary
    return true
end

"""
    periodic_compat(b1::B1, b2::B2)

Checks if two boundaries with the different boundary conditions (B2 and B2) are compatible as opposites. If either is periodic they are not compatible, otherwise there are no restrictions on "matching" boundaries.  
"""
function periodic_compat(b1::B1, b2::B2) where {B1<:Boundary, B2<:Boundary}
    return !(b1.bc isa PeriodicBC || b2.bc isa PeriodicBC)
end

"""
    RectangleDomain{FT}<:AbstractDomain{FT}

Rectangular subtype of AbstractDomain that holds 4 Boundaries and coordinates to form a polygon out of the boundaries. The coordinates field is a RingVec for easy conversion to a linear ring using LibGEOS so that is can be the interior hole of a polygon around the domain. 

In order to create a RectangleDomain, three conditions need to be met. First, if needs to be periodically compatible. This means that pairs of opposite boundaries both need to be periodic if one of them is periodic. Next, the value in the north boundary must be greater than the south boundary and the value in the east boundary must be greater than the west in order to form a valid rectangle. Finally, the last point in coords should match the first one for ease of making into a polygon.
"""
struct RectangleDomain{FT, NBC<:AbstractBC, SBC<:AbstractBC, EBC<:AbstractBC, WBC<:AbstractBC}<:AbstractDomain{FT}
    north::Boundary{NBC, FT}
    south::Boundary{SBC, FT}
    east::Boundary{EBC, FT}
    west::Boundary{WBC, FT}

    RectangleDomain(north, south, east, west) = 
        (periodic_compat(north, south) && periodic_compat(east, west)) &&
        (north.val > south.val && east.val > west.val) ?
        new{typeof(north.val),typeof(north.bc),typeof(south.bc), typeof(east.bc),typeof(west.bc)}(north, south, east, west) : 
        throw(ArgumentError("Periodic boundary must have matching opposite boundary and/or North value must be greater then South and East must be greater than West."))
end

"""
    RectangleDomain(grid, northBC = Open(), southBC = Open(),
    eastBC = Open(), westBC = Open())

Create a maximum rectangular domain from a grid and 4 boundary conditions 
Inputs:
        grid      <Grid>
        southBC   <AbstractBC subtype>
        eastBC    <AbstractBC subtype>
        westBC    <AbstractBC subtype>
        westBC    <AbstractBC subtype>
Output:
        Rectangular Domain with boundaries along the edges of the grid with each boundary having given boundary conditions. 
"""
function RectangleDomain(grid::Grid; northBC = OpenBC(), southBC = OpenBC(),
eastBC = OpenBC(), westBC = OpenBC())
    north = Boundary(northBC, grid.yg[end])
    south = Boundary(southBC, grid.yg[1])
    east = Boundary(eastBC, grid.xg[end])
    west = Boundary(westBC, grid.xg[1])
    return RectangleDomain(north, south, east, west)
end

"""
    CircleDomain{FT}<:AbstractDomain{FT}

Circle subtype of AbstractDomain holds one boundary condition and a Float centroid and radius. The entire domain boundary has the given boundary condition.
"""
struct CircleDomain{FT, BC<:AbstractBC}<:AbstractDomain{FT}
    radius::FT
    centroid::Vector{FT}
    bc::BC

    CircleDomain(radius, centroid, bc) = radius > 0.0 ?
        new{typeof(radius), typeof(bc)}(radius, centroid, bc) :
        throw(ArgumentError("Circle domain radius should be positive."))
end

"""
    CircleDomain(grid::Grid, bc = Open())

Creates a maximum circular domain from a grid and a boundary condition
Inputs:
        grid      <Grid>
        bc        <AbstractBC subtype> 
Output:
        Circular Domain that has a centroid in the center of the grid and a radius of half of the shorter side of the grid. Whole domain has the given boundary condition. 
"""
function CircleDomain(grid::Grid; bc = OpenBC())
    xmin = grid.xg[1]
    ymin = grid.yg[1]
    xrad = (grid.xg[end] - xmin)/2
    yrad = (grid.yg[end] - ymin)/2
    radius = min(xrad, yrad)
    centroid = [xmin + xrad, ymin + yrad]
    return CircleDomain(radius, centroid, bc)
end

"""
    Grid(domain::CircleDomain, nx::Int, ny::Int, t::Type{T} = Float64)

Minimum grid that contains the given CircleDomain with the given number of grid cells in the x and y direction.

Inputs: 
        domain   <CircleDomain>
        dims     <(Int, Int)> grid dimensions - rows -> ny, cols -> nx
        t        <Type> datatype to convert grid fields - must be a Float!
Outputs:
        Square grid where the height and width are equal to domain diameter and centered on circle centroid, with nx grid cells in the x-direction and ny grid cells in the y-direction.
"""
function Grid(domain::CircleDomain, dims::Tuple{Int, Int}, t::Type{T} = Float64) where T
    northval = domain.centroid[2] + domain.radius
    southval = domain.centroid[2] - domain.radius
    eastval = domain.centroid[1] + domain.radius
    westval = domain.centroid[1] - domain.radius
    return Grid(westval, eastval, southval, northval, dims, T)
end

"""
    Grid(domain::RectangleDomain, nx::Int, ny::Int, t::Type{T} = Float64)

Minimum grid that contains the given RecangleDomain with the given number of grid cells in the x and y direction.

Inputs: 
        domain   <RectangleDomain>
        dims     <(Int, Int)> grid dimensions - rows -> ny, cols -> nx
        t        <Type> datatype to convert grid fields - must be a Float!
Outputs:
        Rectangle grid that is exactly the side of the domain with nx grid cells in the x-direction and ny grid cells in the y-direction.
"""
function Grid(domain::RectangleDomain, dims::Tuple{Int, Int},  t::Type{T} = Float64) where T
    northval = domain.north.val
    southval = domain.south.val
    eastval = domain.east.val
    westval = domain.west.val
    return Grid(westval, eastval, southval, northval, dims, T)
end

"""
    Topography{FT<:AbstractFloat}

Singular topographic element with fields describing current state. These are used to create the desired topography within the simulation and will be treated as "islands" within the model in that they will not move or break due to floe interactions, but they will affect floes. 

Coordinates are vector of vector of vector of points of the form:
[[[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]], 
 [[w1, z1], [w2, z2], ..., [wn, zn], [w1, z1]], ...] where the xy coordinates are the exterior border of the element and the wz coordinates, or any other following sets of coordinates, describe holes within it - although there should not be any. This format makes for easy conversion to and from LibGEOS Polygons. 
"""
struct Topography{FT<:AbstractFloat}
    centroid::Vector{FT}
    coords::PolyVec{FT}     # coordinates of topographical element
    height::FT              # height (m)
    area::FT                # area (m^2)
    rmax::FT                # distance of vertix farthest from centroid (m)

    Topography(centroid, coords, height, area, rmax) = 
        height > 0 && area > 0 && rmax > 0 ?
        new{typeof(height)}(centroid, coords, height, area, rmax) :
        throw(ArgumentError("Height, area, and maximum radius of a given topography element should be positive."))
end

"""
    Topography(poly::LG.Polygon, h, t::Type{T} = Float64)

Constructor for topographic element with LibGEOS Polygon
    Inputs:
            poly    <LibGEOS.Polygon> 
            h       <Real> height of element
            t       <Type> datatype to convert ocean fields - must be a Float!
    Output:
            Topographic element with needed fields defined
    Note:
            Types are specified at Float64 below as type annotations given that when written LibGEOS could exclusivley use Float64 (as of 09/29/22). When this is fixed, this annotation will need to be updated.
            We should only run the model with Float64 right now or else we will be converting the Polygon back and forth all of the time. 
"""
function Topography(poly::LG.Polygon, h, t::Type{T} = Float64) where T
    topo = rmholes(poly)
    centroid = LG.GeoInterface.coordinates(LG.centroid(topo))::Vector{Float64}
    area = LG.area(topo)::Float64 
    coords = LG.GeoInterface.coordinates(topo)::PolyVec{Float64}
    rmax = sqrt(maximum([sum(c.^2) for c in translate(coords, -centroid)[1]]))
    return Topography(convert(Vector{T}, centroid), convert(PolyVec{T}, coords),
                      convert(T, h), convert(T, area), convert(T, rmax))
end

"""
    Topography(coords::PolyVec{T}, h, t::Type{T} = Float64)

Topogrpahic element constructor with PolyVec{Float64}(i.e. Vector{Vector{Vector{Float64}}}) coordinates
Inputs:
poly    <LibGEOS.Polygon> 
h       <Real> height of element
t       <Type> datatype to convert ocean fields - must be a Float!
Output:
        Topographic element with needed fields defined
"""
function Topography(coords::PolyVec{<:Real}, h, t::Type{T} = Float64) where T
    return Topography(LG.Polygon(convert(PolyVec{Float64}, coords)), h, T)
    # Polygon convert is needed since LibGEOS only takes Float64 - when this is fixed convert can be removed
end

"""
Singular sea ice floe with fields describing current state. Centroid is a vector of points of the form: [x,y].
Coordinates are vector of vector of vector of points of the form:
[[[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]], 
 [[w1, z1], [w2, z2], ..., [wn, zn], [w1, z1]], ...] where the xy coordinates are the exterior border of the floe
and the wz coordinates, or any other following sets of coordinates, describe holes within the floe.
There should not be holes for the majority of the time as they will be removed, but this format makes for easy
conversion to and from LibGEOS Polygons. 
"""
@kwdef mutable struct Floe{FT<:AbstractFloat}
    centroid::Vector{FT}    # center of mass of floe (might not be in floe!)
    height::FT              # floe height (m)
    area::FT                # floe area (m^2)
    mass::FT                # floe mass (kg)
    moment::FT              # mass moment of intertia
    #angles::Vector{T}
    rmax::FT                # distance of vertix farthest from centroid (m)
    coords::PolyVec{FT}     # floe coordinates
    Δα::FT = 0.0            # change in rotation α since last timestep
    u::FT = 0.0             # floe x-velocity
    v::FT = 0.0             # floe y-velocity
    ξ::FT = 0.0             # floe angular velocity
    fxOA::FT = 0.0          # force from ocean and wind in x direction
    fyOA::FT = 0.0          # force from ocean and wind in y direction
    torqueOA::FT = 0.0      # torque from ocean and wind
    p_dxdt::FT = 0.0        # previous timestep x-velocity
    p_dydt::FT = 0.0        # previous timestep y-velocity
    p_dudt::FT = 0.0        # previous timestep x-acceleration
    p_dvdt::FT = 0.0        # previous timestep x-acceleration
    p_dξdt::FT = 0.0        # previous timestep time derivative of ξ
    p_dαdt::FT = 0.0        # previous timestep ξ
    hflx::FT = 0.0          # heat flux under the floe
    alive::Bool = true      # floe is still active in simulation
end # TODO: do we want to do any checks? Ask Mukund!

"""
    Floe(poly::LG.Polygon, hmean, Δh, ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, t::Type{T} = Float64)

Constructor for floe with LibGEOS Polygon
Inputs:
        poly        <LibGEOS.Polygon> 
        h_mean      <Real> mean height for floes
        Δh          <Real> variability in height for floes
        ρi          <Real> ice density kg/m3 - default 920
        u           <Real> x-velocity of the floe - default 0.0
        v           <Real> y-velcoity of the floe - default 0.0
        ksi         <Real> angular velocity of the floe - default 0.0
        t           <Type> datatype to convert ocean fields - must be a Float!
Output:
        Floe with needed fields defined - all default field values used so all forcings start at 0 and floe is "alive".
        Velocities and the density of ice can be optionally set.
Note:
        Types are specified at Float64 below as type annotations given that when written LibGEOS could exclusivley use Float64 (as of 09/29/22).
        When this is fixed, this annotation will need to be updated.
        We should only run the model with Float64 right now or else we will be converting the Polygon back and forth all of the time. 
"""
function Floe(poly::LG.Polygon, hmean, Δh; ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, t::Type{T} = Float64) where T
    floe = rmholes(poly)
    centroid = LG.GeoInterface.coordinates(LG.centroid(floe))::Vector{Float64}
    h = hmean + (-1)^rand(0:1) * rand() * Δh  # floe height
    area = LG.area(floe)::Float64  # floe area
    mass = area * h * ρi  # floe mass
    moment = calc_moment_inertia(floe, h, ρi = ρi)
    coords = LG.GeoInterface.coordinates(floe)::PolyVec{Float64}
    rmax = sqrt(maximum([sum(c.^2) for c in translate(coords, -centroid)[1]]))

    return Floe(centroid = convert(Vector{T}, centroid),
                height = convert(T, h), area = convert(T, area),
                mass = convert(T, mass), moment = convert(T, moment),
                coords = convert(PolyVec{T}, coords), rmax = convert(T, rmax), u = convert(T, u), v = convert(T, v), ξ = convert(T, ξ))
end

"""
    Floe(coords::PolyVec{Float64}, h_mean, Δh, ρi = 920.0, u = 0.0,
    v = 0.0, ξ = 0.0, t::Type{T} = Float64) where T

Floe constructor with PolyVec{Float64}(i.e. Vector{Vector{Vector{Float64}}}) coordinates
Inputs:
        coords      <Vector{Vector{Vector{Float64}}}> floe coordinates
        h_mean      <Real> mean height for floes
        h_delta     <Real> variability in height for floes
        rho_ice     <Real> ice density kg/m3
        u           <Real> x-velocity of the floe - default 0.0
        v           <Real> y-velcoity of the floe - default 0.0
        ksi         <Real> angular velocity of the floe - default 0.0
        t           <Type> datatype to convert ocean fields - must be a Float!
Output:
        Floe with needed fields defined - all default field values used so all forcings and velocities start at 0 and floe is "alive"
"""
Floe(coords::PolyVec{<:Real}, h_mean, Δh; ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, t::Type{T} = Float64) where T =
    Floe(LG.Polygon(convert(PolyVec{Float64}, coords)), h_mean, Δh, ρi, u, v, ξ, T) 
    # Polygon convert is needed since LibGEOS only takes Float64 - when this is fixed convert can be removed

"""
    domain_in_grid(domain::RectangleDomain, grid)

Checks if given rectangular domain is within given grid and gives user a warning if domain is not of maximum possible size given grid dimensions.

Inputs:
        domain      <RectangularDomain>
        grid        <Grid>
Outputs:
        <Boolean> true if domain is within grid bounds, else false
"""
function domain_in_grid(domain::RectangleDomain, grid)
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
    domain_in_grid(domain::CircleDomain, grid)

Checks if given circular domain is within given grid and gives user a warning if domain is not of maximum possible size given grid dimensions.

Inputs:
        domain      <CircularDomain>
        grid        <Grid>
Outputs:
        <Boolean> true if domain is within grid bounds, else false
"""
function domain_in_grid(domain::CircleDomain, grid)
    x, y = domain.centroid
    r = domain.radius
    gridheight = grid.yg[end] - grid.yg[1]
    gridwidth = grid.xg[end] - grid.xg[1]
    if (x + r <=  grid.xg[end] &&
        x - r >=  grid.xg[1] &&
        y + r <=  grid.yg[end] &&
        y - r >=  grid.yg[1])
        if 2r != min(gridheight, gridwidth)
            @warn "Circle domain does not have maximum possible radius within grid. This could lead to unneeded computation. Consider making the grid smaller of the domain larger."
        end
        return true
    end
    return false
end

@kwdef struct Constants{FT<:AbstractFloat}
    ρi::FT = 920.0              # ice density
    ρo::FT = 1027.0             # ocean density
    ρa::FT = 1.2                # air density
    Cd_io::FT = 3e-3            # ice-ocean drag coefficent
    Cd_ia::FT = 1e-3            # ice-atmosphere drag coefficent
    f::FT = 1.4e-4              # ocean coriolis parameter
    turnθ::FT = 15*pi/180       # ocean turn angle
    L::FT = 2.93e5              # latent heat of freezing [Joules/kg]
    k::FT = 2.14                # Thermal conductivity of surface ice[W/(m*K)]
    A::FT = 70.0                # Upward flux constant A (W/m2)
    B::FT = 10.0                # Upward flux constant B (W/m2/K)
    Q::FT = 200.0               # Solar constant (W/m2)
end

"""
Model which holds grid, ocean, wind structs, each with the same underlying float type (either Float32 of Float64). It also holds an StructArray of floe structs, again each relying on the same underlying float type. Finally it holds several physical constants. These are:
- hflx: difference in ocean and atmosphere temperatures
- h_new: the height of new ice that forms during the simulation
- modulus: elastic modulus used in floe interaction calculations
- ρi: density of ice
- ρo: density of ocean water
- ρa: density of atmosphere
- Cd_io: ice-ocean drag coefficent
- Cd_ia: ice-atmosphere drag coefficent
- f: ocean coriolis forcings
- turn angle: Ekman spiral caused angle between the stress and surface current
              (angle is positive)
"""
struct Model{FT<:AbstractFloat, DT<:AbstractDomain{FT}}
    grid::Grid{FT}
    ocean::Ocean{FT}
    wind::Wind{FT}
    domain::DT
    topos::StructArray{Topography{FT}}
    floes::StructArray{Floe{FT}}
    consts::Constants{FT}
    hflx::Matrix{FT}            # ocean heat flux
    modulus::FT                 # elastic modulus

    Model(grid, ocean, wind, domain, topos, floes, consts, hflx, modulus) =
        (grid.dims == size(ocean.u) == size(wind.u) && domain_in_grid(domain, grid)) ?
        new{typeof(modulus), typeof(domain)}(grid, ocean, wind, domain, topos, floes, consts, hflx, modulus) :
        throw(ArgumentError("Size of grid does not match size of ocean and/or wind OR domain is not within grid."))
end

"""
    Model(grid, ocean, wind, domain, topos, floes, Δt::Int, newfloe_Δt::Int;
    ρi = 920.0, t::Type{T} = Float64)

Model constructor
Inputs:
        grid        <Grid>
        ocean       <Ocean>
        wind        <Wind>
        domain      <AbstractDomain subtype>
        topo        <StructArray{<:Topography}>
        floes       <StructArray{<:Floe}>
        consts      <Contants>
        t           <Type> datatype to convert ocean fields - must
                           be a Float! 
Outputs:
        Model with all needed fields defined and converted to type t.        
"""
function Model(grid, ocean, wind, domain, topos, floes, consts, t::Type{T} = Float64) where T
    hflx = consts.k/(consts.ρi*consts.L) .* (wind.temp .- ocean.temp)
    modulus = 1.5e3*(mean(sqrt.(floes.area)) + minimum(sqrt.(floes.area)))
    return Model(grid, ocean, wind, domain, topos, floes, consts,
                 convert(Matrix{T}, hflx), convert(T, modulus))
end
