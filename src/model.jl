"""
Structs and functions used to define a Subzero model
"""

"""
Grid splitting the model into distinct rectanglular grid cells where xlines are the grid lines in the x-direction (1xn vector) and ylines are the grid lines in the y-direction (1xm vector). xgrid is the xline vector repeated m times as rows in a nxm array and ygrid is the yline vector repeated n times as columns in a nxm vector.
""" #VT<:AbstractVector{T}, MT<:AbstractMatrix{T}
struct Grid{FT<:AbstractFloat}
    size::Tuple{Int, Int}
    xlines::Vector{FT}
    ylines::Vector{FT}
    xgrid::Matrix{FT}
    ygrid::Matrix{FT}
end

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
Singular sea ice floe with fields describing current state. Centroid is a vector of points of the form: [x,y].
Coordinates are vector of vector of vector of points of the form:
[[[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]], 
 [[w1, z1], [w2, z2], ..., [wn, zn], [w1, z1]], ...] where the xy coordinates are the exterior border of the floe and the wz coordinates, or any other following sets of coordinates, describe holes within the floe. There should not be holes for the majority of the time as they will be removed, but this format makes for easy conversion to and from LibGEOS Polygons. 
"""
@kwdef mutable struct Floe{FT<:AbstractFloat}
    centroid::Vector{FT}    # center of mass of floe (might not be in floe!)
    coords::PolyVec{FT}     # floe coordinates
    height::FT              # floe height (m)
    area::FT                # floe area (m^2)
    mass::FT                # floe mass (kg)
    #inertia_moment::T
    #angles::Vector{T}
   
    rmax::FT                # distance of vertix farthest from centroid (m)
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
    coords = translate(LG.GeoInterface.coordinates(floe)::PolyVec{Float64},
                       -centroid)
    αcoords = deepcopy(coords)
    rmax = sqrt(maximum([sum(c.^2) for c in coords[1]]))

    return Floe(centroid = convert(Vector{T}, centroid),
                height = convert(T, h), area = convert(T, a),
                mass = convert(T, m), coords = convert(PolyVec{T}, coords), rmax = convert(T, rmax), αcoords = convert(PolyVec{T}, αcoords),
                ufloe = convert(T, u), vfloe = convert(T, v),
                ξfloe = convert(T, ξ))
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


#TODO add boundary struct
abstract type AbstractBC end

struct Freeflow<:AbstractBC end
struct Periodic<:AbstractBC end  # Not implemented yet
struct Collision<:AbstractBC end # Not implemented yet
struct Compression<:AbstractBC end # Not implemented yet

struct Wall{BC<:AbstractBC, FT<:AbstractFloat}
    condition::BC
    val::FT # this could just be maximum value in x or y direction... 
end

abstract type AbstractBoundary{FT <: AbstractFloat} end

function periodicCompat(w1::W, w2::W) where W<:Wall
    return true
end
function periodicCompat(w1::W1, w2::W2) where {W1<:Wall, W2<:Wall}
    return !(w1.condition isa Periodic || w2.condition isa Periodic)
end


struct RectangleBoundary{FT}<:AbstractBoundary{FT}
    north::Wall
    south::Wall
    east::Wall
    west::Wall
    coords::PolyVec{FT}
    # should we add a coords one that holds the box so we can do polygon operations?
    RectangleBoundary(north, south, east, west, coords) = 
        (periodicCompat(north, south) && periodicCompat(east, west)) &&
        (north.val>south.val && east.val > west.wal) ?
        new{FT}(north, south, east, west, coords) : 
        throw(ArgumentError("Periodic wall must have matching opposite wall and/or North value must be greater then South and East must be greater than West."))
end

function RectangleBoundary(north, south, east, west)
    coords = [[[west.val, north.val], [west.val, south.val],
               [east.val, south.val], [east.val, north.val],
               [west.val, north.val]]]
    return RectangleBoundary(north, south, east, west, coords)
end

# Not implemented yet
struct CircleBoundary{FT, BC<:AbstractBC}<:AbstractBoundary{FT}
    Radius::FT
    Centroid::FT
    condition::BC
    # Should we allow periodic?
end

"""
Model which holds grid, ocean, wind structs, each with the same underlying float type (either Float32 of Float64). It also holds an StructArray of floe structs, again each relying on the same underlying float type. Finally it holds several physical constants. These are:
- heatflux ?
- new_iceh: the height of new ice that forms during the simulation
- ρi: density of ice
- coriolis: ocean coriolis forcings
- turn angle: Ekman spiral caused angle between the stress and surface current
              (angle is positive)
"""
struct Model{FT<:AbstractFloat}
    grid::Grid{FT}
    ocean::Ocean{FT}
    wind::Wind{FT}
    # TODO: Add boundaries
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
function Model(grid, ocean, wind, floes, To, Ta, Δt::Int, newfloe_Δt::Int;
ρi = 920.0, coriolis = 1.4e-4, turnangle = 15*pi/180, L = 2.93e5, k = 2.14, t::Type{T} = Float64) where T
#TODO: Should To and Ta be matricies? They are temperature of Ocean/Ice interface and Atmosphere/Ice interface
    new_iceh, heatflux = calc_new_iceh(To, Ta, Δt, newfloe_Δt, ρi = ρi)
    return Model(grid, ocean, wind, floes, convert(T, heatflux),
                 convert(T, new_iceh), convert(T, ρi),
                 convert(T, coriolis), convert(T, turnangle))
end