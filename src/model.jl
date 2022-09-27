"""
Structs and functions used to define a Subzero model
"""

struct Grid{T<:Real}
    xlines::Vector{T}
    ylines::Vector{T}
    xgrid::Array{T,2}
    ygrid::Array{T,2}
end

"""
if dgrid doesn't evenly divide 2x or 2y you won't get full size
"""
function Grid(x, y, dgrid)
    xlines = collect(-x:dgrid:x) 
    ylines = collect(-y:dgrid:y)
    xgrid = repeat(reshape(xlines, 1, :), inner=(length(ylines),1))
    ygrid = repeat(ylines, outer = (1,length(xlines)))
    return Grid(xlines, ylines, xgrid, ygrid)
end

struct Ocean{T<:Real}
    uocn::Array{T, 2}
    vocn::Array{T, 2}
end

Ocean(grid::Grid) = Ocean(zeros(size(grid.xgrid)), zeros(size(grid.ygrid)))

function Ocean(grid::Grid, psi_func)
    xocn = grid.xgrid
    yocn = grid.ygrid
    dgrid = xocn[1,2] - xocn[1,1]
    psi_ocn = psi_func(xocn, yocn)
    uocn = zeros(size(xocn))
    vocn = zeros(size(yocn))
    # TODO: I don't know if this is standard, ask Mukund
    uocn[2:end,:] = -(psi_ocn[2:end,:]-psi_ocn[1:end-1,:])/dgrid
    vocn[:,2:end]=(psi_ocn[:,2:end]-psi_ocn[:,1:end-1])/dgrid
    return Ocean( uocn, vocn)
end

struct Wind{T<:Real}
    uwind::Array{T, 2}
    vwind::Array{T, 2}
end

Wind(grid::Grid) =  Wind(zeros(size(grid.xgrid)), zeros(size(grid.ygrid)))

Wind(ocn::Ocean) = Wind(ocn.uocn, ocn.vocn)



mutable struct Floe{T<:Real}
    xfloe::T  # centroid x-coordinate
    yfloe::T  # centroid y-coordinate
    hfloe::T  # floe height
    ufloe::T  # floe x-velocity
    vfloe::T  # floe y-velocity
    ksifloe::T  # floe angular velocity
    coords::Vector{Vector{Vector{T}}}  # floe coordinates
    alpha::T  # angle rotated from reference coordinates coords
    coords_alpha::Vector{Vector{Vector{T}}}  # rotated coordinates
    inertia_moment::T  
    area::T
    mass::T
    rmax::T  # distance of vertix farthest from the centroid (m)
    fxOA::T  # force from ocean and atmosphere (wind) in x direction on floe
    fyOA::T  # force from ocean and atmosphere (wind) in y direction on floe
    torqueOA::T  # torque from ocean and atmosphere (wind) on floe
end

struct Model
    ocean::Ocean
    wind::Wind
    floes::StructArray{Floe}
end


