# # Basic plotting functions (and stub functions) for Subzero simulations
export plot_sim, plot_sim_with_ocean_field, prettytime, get_curl, calc_ro_field

#= 
## What plotting functionality is availible for Subzero simulations?

We provide very basic plotting functionalities for Subzero simulations. Given plotting is so
customizable, we just wanted to provide a roadmap for how a user might plot a simulation and
provide a few utility functions.

Below are the stubs of functions that require `CairoMakie`. The actual code for these
functions is in `ext/SubzeroMakieExt`. These functions are only loaded when the user loads
`CairoMakie`. This prevents unneccesarily long compiling times when the user just wants to
run simulations. They can then load `CairoMakie` once post simulation runs and create all of
the needed videos at once, rather than loading `CairoMakie` prior to each simulation.

We also provide a few basic plotting utility functions, like `prettytime` here as they do
not depend on `CairoMakie` and therefore don't need to be abstracted away into the
extension.
=#

# Stub functions that depend on CairoMakie implemented in ext/SubzeroMakieExt.jl

function plot_sim end
function plot_sim_with_ocean_field end

# Constants used in plotting code
const FLOE_FN_DEF = "floe outputwriter output file path and name"
const INITIAL_STATE_FN_DEF = "initial state outputwriter output file path and name"
const OCEAN_FN_DEF = "`Oceananigans` output surface.nc file"
const ΔT_DEF = "length of timestep in integer seconds"
const MP4_OUTPUT_FN = "output video file path and name (should end with .mp4)"

# Utility functions (and example functions) that don't depend on CairoMakie

"""
    prettytime(t)

Turn time in seconds into units of minutes, hours, days, or years as appropriate
## Aguments:
- `t::Real`: number of seconds
## Returns:
- `::String`: number of seconds converted to a string value in minutes, hours, days, or years with units
## Note:
    This code was modified from the this [source code](https://github.com/JuliaCI/BenchmarkTools.jl/blob/master/src/trials.jl).
"""
function prettytime(t)
    minute = 60
    hour = 60 * minute
    day = 24 * hour
    year = 365 * day

    iszero(t) && return "0 seconds"
    if t < minute
        value = floor(Int, t)
        units = value == 1 ? "second" : "seconds"
    elseif t < hour
        value = floor(Int, t / minute)
        units = value == 1 ? "minute" : "minutes"
    elseif t < day
        value = floor(Int, t / hour)
        units = value == 1 ? "hour" : "hours"
    elseif t < year
        value = floor(Int, t / day)
        units = value == 1 ? "day" : "days"
    else
        value = floor(Int, t / year)
        units = value == 1 ? "year" : "years"
    end
    return @sprintf("%d %s", value, units)
end

"""
    get_curl(fldx,fldy,dx,dy)

Calculate curl from ocean u and v velocity fields
## Arguments:
-`fldx::Matrix{AbstractFloat}`: ocean `u` velocity field
-`fldy::Matrix{AbstractFloat}`: ocean `v` velocity field
-`dx::AbstractFloat`:: x-distance over which velocity fields are provided
-`dy::AbstractFloat`: y-distance over which velocity fields are provided
## Returns:
- `::Matrix{AbstractFloat}`: ocean curl at the center of each grid cell
"""
function get_curl(fldx,fldy,dx,dy)
    # fldx must be on the u-grid point and fldy on the v-grid
    # returns field on the omega grid
    nx,ny = size(fldx)
    dvdx = zeros(ny,nx)
    dudy = zeros(ny,nx)
    dvdx[1:end-1,:] = diff(fldy,dims=1)/dx
    dudy[:,1:end-1] = diff(fldx,dims=2)/dy
    return (dvdx - dudy)
end

"""
    calc_ro_field(ocean_fn)

Calculate surface vorticity from an Oceananigan's ocean file.

## Arguments:
- `ocean_fn::String`: $OCEAN_FN_DEF
Returns:
- `ro::Array{AbstractFloat}`: 3D array where first two dimensions are ocean size (Nx, Ny) and the third dimension is time over the simulation
- `xc::Vector{AbstractFloat}`: x grid points for ro values
- `yc::Vector{AbstractFloat}`: y grid points for ro values
"""
function calc_ro_field(ocean_fn)
    xc = NetCDF.ncread(ocean_fn, "xC")
    yc = NetCDF.ncread(ocean_fn, "yC")
    dx = xc[2] - xc[1]
    dy = yc[2] - yc[1]
    Nx = length(xc)
    Ny = length(yc)
    usurf = NetCDF.ncread(ocean_fn, "u")[1:Nx,1:Ny,:,:]
    vsurf = NetCDF.ncread(ocean_fn, "v")[1:Nx,1:Ny,:,:]
    omega = 2*π/(3600*24)
    f = f = 2*omega*sin(70*π/180)
    nsteps = size(usurf, 4)
    ro = zeros(Nx, Ny, nsteps)
    @views for i in 1:nsteps
        ro[:, :, i] .= get_curl(usurf[:,:,1,i], vsurf[:,:,1,i], dx,dy) ./ f
    end
    return ro, xc, yc
end