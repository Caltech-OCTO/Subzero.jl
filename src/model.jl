"""
Structs and functions used to define a Subzero model
"""

"""
    Grid{FT<:AbstractFloat}

Grid splitting the model into distinct rectanglular grid cells where xlines are the grid lines in the x-direction (1xn vector) and ylines are the grid lines in the y-direction (1xm vector). xgrid is the xline vector repeated m times as rows in a nxm array and ygrid is the yline vector repeated n times as columns in a nxm vector.
"""
struct Grid{FT<:AbstractFloat}
    size::Tuple{Int, Int}
    xlines::Vector{FT}
    ylines::Vector{FT}
    xgrid::Matrix{FT}
    ygrid::Matrix{FT}
end
# TODO: Add check that ocean/wind are the same size as the grid

"""
    Grid(x, y, dgrid)

Construct a rectanglular grid for the model.
Inputs: 
        x       <Real> grid length will range from -x to x
        y       <Real> grid height will range from -y to ygrid
        Δgrid   <Real> length/height of square grid cells
        t       <Type> datatype to convert grid fields - must be a Float!
Output: 
        Grid with length of 2x and height of 2y with square dgrid x dgrid grid cells
Warning: If dgrid doesn't evenly divide 2x or 2y you won't get full size grid.
         The grid will be "trimmed" to the nearest full grid square.
"""
function Grid(x, y, Δgrid, t::Type{T} = Float64) where T
    xlines = collect(T, -x:Δgrid:x) 
    ylines = collect(T, -y:Δgrid:y)
    nxlines = length(xlines)
    nylines = length(ylines)
    xgrid = repeat(reshape(xlines, 1, :), inner=(nylines,1))
    ygrid = repeat(ylines, outer = (1, nxlines))
    return Grid((nxlines, nylines), xlines, ylines, xgrid, ygrid)
end

"""
Ocean velocities in the x-direction (uocn) and y-direction (vocn). uocn and vocn should match the size of the corresponding model grid so that there is one x and y velocity value for each grid cell. Ocean also needs temperature in each grid cell. Model cannot be constructed if size of ocean and grid do not match.
"""
struct Ocean{FT<:AbstractFloat}
    uocn::Matrix{FT}
    vocn::Matrix{FT}
    tempocn::Matrix{FT}
end

"""
    Ocean(grid, u, v, temp, FT)

Construct model ocean.
Inputs: 
        grid    <Grid> model grid cell
        u       <Real> ocean x-velocity for each grid cell
        v       <Real> ocean y-velocity for each grid cell
        temp    <Real> ocean temperature for each grid cell
        t       <Type> datatype to convert ocean fields - must be a Float!
Output: 
        Ocean with constant velocity and temperature in each grid cell.
"""
Ocean(grid, u, v, temp, t::Type{T} = Float64) where T =
    Ocean(fill(convert(T, u), grid.size), 
          fill(convert(T, v), grid.size), 
          fill(convert(T, temp), grid.size))

# TODO: Do we want to be able to use a psi function? - Ask Mukund

"""
Wind velocities in the x-direction (uwind) and y-direction (vwind). uwind and vwind should match the size of the corresponding model grid so that there is one x and y velocity value for each grid cell. Model cannot be constructed if size of wind and grid do not match.
"""
struct Wind{FT<:AbstractFloat}
    uwind::Matrix{FT}
    vwind::Matrix{FT}
end

"""
    Wind(grid, u, v, FT)

Construct model atmosphere/wind.
Inputs: 
        grid    <Grid> model grid cell
        u       <Real> wind x-velocity for each grid cell
        v       <Real> wind y-velocity for each grid cell
        t       <Type> datatype to convert ocean fields - must be a Float!
Output: 
        Ocean with constant velocity and temperature in each grid cell.
"""
Wind(grid, u, v, t::Type{T} = Float64) where T = 
    Wind(fill(convert(T, u), grid.size), fill(convert(T, v), grid.size))

"""
    Wind(ocn)

Construct model atmosphere/wind.
Inputs: 
        ocn     <Ocean> model ocean
Output: 
        Ocean with the same wind velocity as ocean velocity in each grid cell.
Note:
        Type does not need to be set at the ocean type will already be set.
"""
Wind(ocn) = Wind(deepcopy(ocn.uocn), deepcopy(ocn.vocn))


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
    condition::BC
    val::FT # this could just be maximum value in x or y direction... 
end

"""
    AbstractDomain{FT <: AbstractFloat}

An abstract type for shapes of domains that contain the model.
"""
abstract type AbstractDomain{FT <: AbstractFloat} end

"""
    periodicCompat(b1::B, b2::B)

Checks if two boundaries with the same boundary condition B are compatible as opposites. They are by definition, given that if both floes are periodic they work as a north/south pair or a east/west pair and otherwise there are no restrictions on "matching" boundaries.
"""
function periodicCompat(b1::B, b2::B) where B<:Boundary
    return true
end

"""
    periodicCompat(b1::B1, b2::B2)

Checks if two boundaries with the different boundary conditions (B2 and B2) are compatible as opposites. If either is periodic they are not compatible, otherwise there are no restrictions on "matching" boundaries.  
"""
function periodicCompat(b1::B1, b2::B2) where {B1<:Boundary, B2<:Boundary}
    return !(b1.condition isa PeriodicBC || b2.condition isa PeriodicBC)
end

"""
    RectangleDomain{FT}<:AbstractDomain{FT}

Rectangular subtype of AbstractDomain that holds 4 Boundaries and coordinates to form a polygon out of the boundaries. The coordinates field is a RingVec for easy conversion to a linear ring using LibGEOS so that is can be the interior hole of a polygon around the domain. 

In order to create a RectangleDomain, three conditions need to be met. First, if needs to be periodically compatible. This means that pairs of opposite boundaries both need to be periodic if one of them is periodic. Next, the value in the north boundary must be greater than the south boundary and the value in the east boundary must be greater than the west in order to form a valid rectangle. Finally, the last point in coords should match the first one for ease of making into a polygon.
"""
struct RectangleDomain{FT}<:AbstractDomain{FT}
    north::Boundary
    south::Boundary
    east::Boundary
    west::Boundary
    coords::RingVec{FT} 

    RectangleDomain(north, south, east, west, coords) = 
        (periodicCompat(north, south) && periodicCompat(east, west)) &&
        (north.val>south.val && east.val > west.val) ?
        new{eltype(eltype(coords))}(north, south, east, west,
                                    valid_ringvec!(coords)) : 
        throw(ArgumentError("Periodic boundary must have matching opposite boundary and/or North value must be greater then South and East must be greater than West."))
end

"""
    RectangleDomain(north, south, east, west)

Create a rectangular domain from 4 boundaries. 
Inputs:
        north   <Boundary>
        south   <Boundary>
        east    <Boundary>
        west    <Boundary>
Output:
        Rectangular Domain with 4 given boundaries and rectangle coordinates created by using each boundaties values to determine the corners.
"""
function RectangleDomain(north::Boundary, south::Boundary, east::Boundary, west::Boundary)
    coords = [[west.val, north.val], [west.val, south.val],
              [east.val, south.val], [east.val, north.val],
              [west.val, north.val]]
    return RectangleDomain(north, south, east, west, coords)
end

"""
    RectangleDomain(grid, northBC = Open(), southBC = Open(),
    eastBC = Open(), westBC = Open())

Create a rectangular domain from a grid and 4 boundary conditions 
Inputs:
        grid      <Grid>
        southBC   <AbstractBC subtype>
        eastBC    <AbstractBC subtype>
        westBC    <AbstractBC subtype>
        westBC    <AbstractBC subtype>
Output:
        Rectangular Domain with boundaries along the edges of the grid with each boundary having given boundary conditions. 
"""
function RectangleDomain(grid::Grid, northBC = Open(), southBC = Open(),
eastBC = Open(), westBC = Open())
    northwall = Boundary(northBC, maximum(grid.ylines))
    southwall = Boundary(southBC, minimum(grid.ylines))
    eastwall = Boundary(eastBC, maximum(grid.xlines))
    westwall = Boundary(westBC, minimum(grid.xlines))
    return RectangleDomain(northwall, southwall, eastwall, westwall)
end

"""
    CircleDomain{FT}<:AbstractDomain{FT}

Circle subtype of AbstractDomain holds one boundary condition and a Float centroid and radius. The entire domain boundary has the given boundary condition.

Note: Not implemented yet!!
"""
struct CircleDomain{FT, BC<:AbstractBC}<:AbstractDomain{FT}
    Radius::FT
    Centroid::FT
    condition::BC
end

struct Topography{FT<:AbstractFloat}
    centroid::Vector{FT}    # center of mass of topographical element
    coords::PolyVec{FT}     # coordinates
    height::FT              # height (m)
    area::FT                # area (m^2)
    rmax::FT                # distance of vertix farthest from centroid (m)
end

function Topography(poly::LG.Polygon, h, t::Type{T} = Float64) where T
    topo = rmholes(poly)
    centroid = LG.GeoInterface.coordinates(LG.centroid(topo))::Vector{Float64}
    a = LG.area(topo)::Float64  # floe area
    coords = translate(LG.GeoInterface.coordinates(topo)::PolyVec{Float64},
                       -centroid)
    rmax = sqrt(maximum([sum(c.^2) for c in coords[1]]))
    return Topography(centroid, coords, h, a, rmax)
end

function Topography(coords::PolyVec{T}, h, t::Type{T} = Float64) where T
    return Topography(LG.Polygon(coords), h, T)
end

"""
Singular sea ice floe with fields describing current state. Centroid is a vector of points of the form: [x,y].
Coordinates are vector of vector of vector of points of the form:
[[[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]], 
 [[w1, z1], [w2, z2], ..., [wn, zn], [w1, z1]], ...] where the xy coordinates are the exterior border of the floe and the wz coordinates, or any other following sets of coordinates, describe holes within the floe. There should not be holes for the majority of the time as they will be removed, but this format makes for easy conversion to and from LibGEOS Polygons. 
"""
@kwdef mutable struct Floe{FT<:AbstractFloat}
    centroid::Vector{FT}    # center of mass of floe (might not be in floe!)
    height::FT              # floe height (m)
    area::FT                # floe area (m^2)
    mass::FT                # floe mass (kg)
    moment::FT      # mass moment of intertia
    #angles::Vector{T}
    rmax::FT                # distance of vertix farthest from centroid (m)
    coords::PolyVec{FT}     # floe coordinates
    αcoords::PolyVec{FT}    # rotated coordinates (rotated by angle alpha)
    α::FT = 0.0             # angle rotated from starting coords
    ufloe::FT = 0.0         # floe x-velocity
    vfloe::FT = 0.0         # floe y-velocity
    ξfloe::FT = 0.0         # floe angular velocity
    fxOA::FT = 0.0          # force from ocean and wind in x direction
    fyOA::FT = 0.0          # force from ocean and wind in y direction
    torqueOA::FT = 0.0      # torque from ocean and wind
    alive::Bool = true      # floe is still active in simulation
end

"""
    Floe(poly::LG.Polygon, hmean, Δh, ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, t::Type{T} = Float64) where T

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
        Floe with needed fields defined - all default field values used so all forcings start at 0 and floe is "alive". Velocities and the density of ice can be optionally set.
Note:
        Types are specified at Float64 below as type annotations given that when written LibGEOS could exclusivley use Float64 (as of 09/29/22). When this is fixed, this annotation will need to be updated.
        We should only run the model with Float64 right now or else we will be converting the Polygon back and forth all of the time. 
"""
function Floe(poly::LG.Polygon, hmean, Δh, ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, t::Type{T} = Float64) where T
    floe = rmholes(poly)
    centroid = LG.GeoInterface.coordinates(LG.centroid(floe))::Vector{Float64}
    h = hmean + (-1)^rand(0:1) * rand() * Δh  # floe height
    a = LG.area(floe)::Float64  # floe area
    m = a * h * ρi  # floe mass
    mi = calc_moment_inertia(floe, h, ρi = ρi)
    coords = translate(LG.GeoInterface.coordinates(floe)::PolyVec{Float64},
                       -centroid)
    αcoords = deepcopy(coords)
    rmax = sqrt(maximum([sum(c.^2) for c in coords[1]]))

    return Floe(centroid = convert(Vector{T}, centroid),
                height = convert(T, h), area = convert(T, a),
                mass = convert(T, m), moment = convert(T, mi),
                coords = convert(PolyVec{T}, coords), rmax = convert(T, rmax), αcoords = convert(PolyVec{T}, αcoords), ufloe = convert(T, u), vfloe = convert(T, v), ξfloe = convert(T, ξ))
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
Floe(coords::PolyVec{Float64}, h_mean, Δh, ρi = 920.0, u = 0.0,
v = 0.0, ξ = 0.0, t::Type{T} = Float64) where T =
    Floe(LG.Polygon(coords), h_mean, Δh, ρi, u, v, ξ, T)

"""
Model which holds grid, ocean, wind structs, each with the same underlying float type (either Float32 of Float64). It also holds an StructArray of floe structs, again each relying on the same underlying float type. Finally it holds several physical constants. These are:
- heatflux ?
- new_iceh: the height of new ice that forms during the simulation
- ρi: density of ice
- coriolis: ocean coriolis forcings
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
    heatflux::FT                # ocean heat flux
    new_iceh::FT                # thickness of new sea ice during model run
    ρi::FT                      # ice density
    coriolis::FT                # ocean coriolis force
    turnangle::FT               # ocean turn angle
end

"""
    calc_new_iceh(To, Ta, Δt, newfloe_Δt; ρi = 920.0, L = 2.93e5, k = 2.14)

Calculate height of new floes during simulation and heat flux.
Inputs:
        To          <Real> Ocean-Ice interface temperature
        Ta          <Real> Atmosphere-Ice interface temperature
        Δt          <Int> Timestep within the model in seconds
        newfloe_Δt  <Int> Timesteps between new floe creation
        ρi          <Float> Density of ice
        coriolis    <Float> Ocean coriolis forcings
        turnangle   <Float> Ekman spiral caused angle between the stress and
                            surface current - angle is positive
        L           <Float> Latent heat of freezing [Joules/kg]
        k           <Float> Thermal conductivity of surface ice
                            [Watts/(meters Kelvin)]
Outputs:
        heatflux: ?
        h0: the height of new ice that forms during the simulation
"""
function calc_new_iceh(To, Ta, Δt, newfloe_Δt; ρi = 920.0, L = 2.93e5, k = 2.14)
    h0 = real(sqrt(2k * Δt * newfloe_Δt * (To - Ta)/ρi*L))
    heatflux = k./h0.*(To - Ta)
    return h0, heatflux
end

"""
    Model(grid, ocean, wind, floe, To, Ta, Δt::Int, newfloe_Δt::Int;
    ρi = 920.0, coriolis = 1.4e-4, turnangle = 15*pi/180, L = 2.93e5, k = 2.14, t::Type{T} = Float64) where T

Model constructor
Inputs:
        grid        <Grid>
        ocean       <Ocean>
        wind        <Wind>
        domain      <AbstractDomain subtype>
        topo        <StructArray{<:Topography}>
        floes       <StructArray{<:Floe}>
        To          <Real> Ocean-Ice interface temperature
        Ta          <Real> Atmosphere-Ice interface temperature
        Δt          <Int> Timestep within the model in seconds
        newfloe_Δt  <Int> Timesteps between new floe creation
        ρi          <Float> Density of ice
        coriolis    <Float> Ocean coriolis forcings
        turnangle   <Float> Ekman spiral caused angle between the stress and
                            surface current - angle is positive
        L           <Float> Latent heat of freezing [Joules/kg]
        k           <Float> Thermal conductivity of surface ice
                            [Watts/(meters Kelvin)]
        t           <Type> datatype to convert ocean fields - must
                           be a Float! 
Outputs:
        Model with all needed fields defined.         
"""
function Model(grid, ocean, wind, domain, topo, floes, To, Ta, Δt::Int, newfloe_Δt::Int; ρi = 920.0, coriolis = 1.4e-4, turnangle = 15*pi/180,
L = 2.93e5, k = 2.14, t::Type{T} = Float64) where T
#TODO: Should To and Ta be matricies? They are temperature of Ocean/Ice interface and Atmosphere/Ice interface
    new_iceh, heatflux = calc_new_iceh(To, Ta, Δt, newfloe_Δt, ρi = ρi)
    return Model(grid, ocean, wind, domain, topo, floes, convert(T, heatflux),
                 convert(T, new_iceh), convert(T, ρi),
                 convert(T, coriolis), convert(T, turnangle))
end