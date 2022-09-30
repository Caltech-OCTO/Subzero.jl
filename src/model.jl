"""
Structs and functions used to define a Subzero model
"""

"""
Grid splitting the model into distinct rectanglular grid cells where xlines are the grid lines in the x-direction (1xn vector) and ylines are the grid lines in the y-direction (1xm vector). xgrid is the xline vector repeated m times as rows in a nxm array and ygrid is the yline vector repeated n times as columns in a nxm vector.
"""
struct Grid{T<:AbstractFloat, VT<:AbstractVector{T}, MT<:AbstractMatrix{T}}
    size::Tuple{Int, Int}
    xlines::VT
    ylines::VT
    xgrid::MT
    ygrid::MT
end

"""
    Grid(x, y, dgrid)

Construct a rectanglular grid for the model.
Inputs: 
        x       <Real> grid length will range from -x to x
        y       <Real> grid height will range from -y to ygrid
        dgrid   <Real> length/height of square grid cells
        t       <Type> datatype to convert grid fields - must be a Float!
Output: 
        Grid with length of 2x and height of 2y with square dgrid x dgrid grid cells
Warning: If dgrid doesn't evenly divide 2x or 2y you won't get full size grid.
         The grid will be "trimmed" to the nearest full grid square.
"""
function Grid(x, y, dgrid, t::Type{T} = Float64) where T
    xlines = collect(T, -x:dgrid:x) 
    ylines = collect(T, -y:dgrid:y)
    nxlines = length(xlines)
    nylines = length(ylines)
    xgrid = repeat(reshape(xlines, 1, :), inner=(nylines,1))
    ygrid = repeat(ylines, outer = (1, nxlines))
    return Grid((nxlines, nylines), xlines, ylines, xgrid, ygrid)
end

"""
Ocean velocities in the x-direction (uocn) and y-direction (vocn). uocn and vocn should match the size of the corresponding model grid so that there is one x and y velocity value for each grid cell. Ocean also needs temperature in each grid cell. Model cannot be constructed if size of ocean and grid do not match.
"""
struct Ocean{T<:AbstractFloat, MT<:AbstractMatrix{T}}
    uocn::MT
    vocn::MT
    tempocn::MT
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
struct Wind{T<:AbstractFloat, MT<:AbstractMatrix{T}}
    uwind::MT
    vwind::MT
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
@kwdef mutable struct Floe{T<:AbstractFloat, VT<:AbstractVector{T}, PT<:AbstractVector{<:AbstractVector{VT}}}
    centroid::VT      # center of mass of floe (might not be in floe!)
    coords::PT         # floe coordinates
    height::T                   # floe height (m)
    area::T                     # floe area (m^2)
    mass::T                     # floe mass (kg)
    #inertia_moment::T
    #angles::Vector{T}
   
    rmax::T                     # distance of vertix farthest from centroid (m)
    coords_alpha::PT    # rotated coordinates (rotated by angle alpha)
    alpha::T = 0.0              # angle rotated from starting coords
    ufloe::T = 0.0              # floe x-velocity
    vfloe::T = 0.0              # floe y-velocity
    ksifloe::T = 0.0            # floe angular velocity
    fxOA::T = 0.0               # force from ocean and wind in x direction
    fyOA::T = 0.0               # force from ocean and wind in y direction
    torqueOA::T = 0.0           # torque from ocean and wind
    alive::Bool = true          # floe is still active in simulation
end

"""
Constructor for floe with LibGEOS Polygon
Inputs:
        poly        <LibGEOS.Polygon> 
        h_mean      <Real> mean height for floes
        h_delta     <Real> variability in height for floes
        rho_ice     <Real> ice density kg/m3 - default 920
        u           <Real> x-velocity of the floe - default 0.0
        v           <Real> y-velcoity of the floe - default 0.0
        ksi         <Real> angular velocity of the floe - default 0.0
Output:
        Floe with needed fields defined - all default field values used so all forcings start at 0 and floe is "alive". Velocities and the density of ice can be optionally set.
Note:
        Types are specified at Float64 below as type annotations given that when written LibGEOS could exclusivley use Float64 (as of 09/29/22). When this is fixed, this annotation will need to be updated.
        We should only run the model with Float64 right now or else we will be converting the Polygon back and forth all of the time. 
"""
function Floe(poly::LG.Polygon, hmean, hdelta, rho = 920.0, u = 0.0, v = 0.0, ksi = 0.0, t::Type{T} = Float64) where T
    floe = rmholes(poly)
    centroid = LG.GeoInterface.coordinates(LG.centroid(floe))::Vector{Float64}
    height = hmean + (-1)^rand(0:1) * rand() * hdelta
    A = LG.area(floe)::Float64
    mass = A*height*rho
    coords = translate(
            LG.GeoInterface.coordinates(floe)::PolyVec{Float64},
            -centroid)
    coords_alpha = deepcopy(coords)
    rmax = sqrt(maximum([sum(c.^2) for c in coords[1]]))

    return Floe(centroid = convert(Vector{T}, centroid),
                height = convert(T, height), area = convert(T, A),
                mass = convert(T, mass), coords = convert(PolyVec{T}, coords), rmax = convert(T, rmax),
                coords_alpha = convert(PolyVec{T}, coords_alpha),
                ufloe = convert(T, u), vfloe = convert(T, v),
                ksifloe = convert(T, ksi))
end

"""
Floe constructor with PolyVec{Float64}(i.e. Vector{Vector{Vector{Float64}}}) coordinates
Inputs:
        coords      <Vector{Vector{Vector{Float64}}}> floe coordinates
        h_mean      <Real> mean height for floes
        h_delta     <Real> variability in height for floes
        rho_ice     <Real> ice density kg/m3
        u           <Real> x-velocity of the floe - default 0.0
        v           <Real> y-velcoity of the floe - default 0.0
        ksi         <Real> angular velocity of the floe - default 0.0
Output:
        Floe with needed fields defined - all default field values used so all forcings and velocities start at 0 and floe is "alive"
"""
Floe(coords::PolyVec{Float64}, h_mean, h_delta, rho_ice = 920.0, u = 0.0,
v = 0.0, ksi = 0.0, t::Type{T} = Float64) where T =
    Floe(LG.Polygon(coords), h_mean, h_delta, rho_ice, u, v, ksi, T)


# TODO: Add boundary struct, constructors, update Model struct to take in correct specific types, update model constructor to convert to type T

struct Model{T<:AbstractFloat, VT<:AbstractVector{T},
MT<:AbstractMatrix{T}, PT<:AbstractVector{<:AbstractVector{VT}},
GT<:Grid{T, VT, MT}, OT<:Ocean{T, MT}, WT<:Wind{T, MT}, FT<:Floe{T, VT, PT}}
# TODO: Add boundaries
    grid::GT
    ocean::OT
    wind::WT
    floe::FT
    mean_heatflux::T           # ocean heat flux
    iceheight::T               # thickness of new sea ice during model run
    rhoice::T                  # ice density
    coriolis::T                # ocean coriolis force
    turnangle::T               # ocean turn angle
end

function calc_icestuff(ocean, icetemp, dt, newfloe_dt, t::Type{T} = Float64)where T
    heatflux = (7.4*10^(-8)/72).*(icetemp .- ocean.tempocn)
    mean_heatflux = mean(heatflux)
    iceheight = -2*dt*mean_heatflux*newfloe_dt
    return convert(T, iceheight), convert(T, mean_heatflux)
end

function Model(grid, ocean, wind, floe, icetemp, dt::Int, newfloe_dt::Int,
               rhoice = 920.0, coriolis = 1.4e-4, turnangle = 15*pi/180,
               t::Type{T} = Float64)where T

    heatflux = (7.4*10^(-8)/72).*(icetemp .- ocean.tempocn)
    mean_heatflux = mean(heatflux)
    iceheight = -2*dt*mean_heatflux*newfloe_dt
    return Model(grid, ocean, wind, floe, convert(T, mean_heatflux),
                 convert(T, iceheight), convert(T, rhoice),
                 convert(T, coriolis), convert(T, turnangle))
end
