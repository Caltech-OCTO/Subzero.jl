"""
Structs and functions used to define a Subzero oceans
"""

"""
    IceStressCell{FT<:AbstractFloat}

Struct to collect stress from ice floes on ocean grid cells. One IceStressCell
corresponds to one grid cell. It holds a list of running totals of stress on the
cell, and a running list of the number of points making up those running totals.
Each element in the list corresponds to one floe, which is denoted in the
corresponding CellFloes matrix within the grid. 
"""
struct IceStressCell{FT<:AbstractFloat}
    τx::Vector{FT}
    τy::Vector{FT}
    npoints::Vector{Int}
end

"""
    IceStressCell{FT}()

Constructs an IceStressCell object with empty lists for fields.
"""
IceStressCell{FT}() where {FT} = IceStressCell{FT}(
    Vector{FT}(),
    Vector{FT}(),
    Vector{Int}()
)
"""
    empty!(cell::IceStressCell)

Empties the vectors within an IceStressCell
"""
function Base.empty!(cell::IceStressCell)
    empty!(cell.τx)
    empty!(cell.τy)
    empty!(cell.npoints)
    return
end

"""
    Ocean{FT<:AbstractFloat}

Simulation ocean holding ocean values on the grid-scale with matricies of the
same size as the model's grid. The struct has the following fields:
- u is the ocean velocities in the x-direction for each grid cell
- v is the ocean velocities in the y-direction for each grid cell
- temp is the ocean temperature for each grid cell
- hflx_factor is a factor to calculate the ocean-atmosphere heat flux for a 
  cell in that grid cell by multiplying by its height
- si_frac is the fraction of area in each grid cell that is covered in sea-ice

Ocean fields must all be matricies with dimensions equal to the number of grid
lines in the model's grid. 
Note: If a periodic boundary is used in the domain, the last grid cell in that
direction will not be used as it is equivalent to the first grid cell. Thus, for
all of these fields, the first and last value in the x and/or y direction should
be equal if the east-west or north-south boundary pair are periodic
respectively.
"""
struct Ocean{FT<:AbstractFloat}
    u::Matrix{FT}
    v::Matrix{FT}
    temp::Matrix{FT}
    hflx_factor::Matrix{FT}
    scells::Matrix{IceStressCell{FT}}
    τx::Matrix{FT}
    τy::Matrix{FT}
    si_frac::Matrix{FT}
    dissolved::Matrix{FT}

    function Ocean{FT}(
        u,
        v,
        temp,
        hflx,
        scells,
        τx,
        τy,
        si_frac,
        dissolved,
    ) where {FT <: AbstractFloat}
        if !all(-3 .<= temp .<= 0)
            @warn "Ocean temperatures are above the range for freezing. The \
                thermodynamics aren't currently setup for these conditions."
        end
        if !(size(u) == size(v) == size(τx) == size(τy))
            throw(ArgumentError("One or more of the ocean vector fields aren't \
                the same dimension."))
        end
        if !(size(temp) == size(hflx) == size(si_frac) == size(dissolved))
            throw(ArgumentError("One or more of the ocean tracer fields aren't \
                the same dimension."))
        end
        new{FT}(u, v, temp, hflx, scells, τx, τy, si_frac, dissolved)
    end
end

"""
    Ocean(::Type{FT}, args...)

A float type FT can be provided as the first argument of any Ocean constructor.
An Ocean of type FT will be created by passing all other arguments to the
correct constructor. 
"""
Ocean(::Type{FT}, args...) where {FT <: AbstractFloat} =
    Ocean{FT}(args...)

"""
    Ocean(args...)

If a type isn't specified, Ocean will be of type Float64 and the correct
constructor will be called with all other arguments.
"""
Ocean(args...) = Ocean{Float64}(args...)

"""
    Ocean{FT}(u, v, temp)

Construct model ocean.
Inputs:
    u       <Matrix> ocean x-velocity matrix with u for each grid line
    v       <Matrix> ocean y-velocity matrix with u for each grid line
    temp    <Matrix> temperature matrix with ocean/ice interface temperature for
                each grid line
Output: 
    Model ocean with given velocity and temperature fields on each grid line.
"""
function Ocean{FT}(
    u,
    v,
    temp,
) where {FT <: AbstractFloat}
    Nx, Ny = size(temp) .- 1
    return Ocean{FT}(
        u,
        v,
        temp,
        zeros(FT, Nx + 1, Ny + 1), # heat flux
        [IceStressCell{FT}() for i in 1:(Nx + 1), j in 1:(Ny + 1)],
        zeros(FT, Nx + 1, Ny + 1),  # x-stress
        zeros(FT, Nx + 1, Ny + 1),  # y-stress
        zeros(FT, Nx + 1, Ny + 1),  # sea ice fraction
        zeros(FT, Nx + 1, Ny + 1),  # dissolved
    )
end

"""
    Ocean{FT}(grid, u, v, temp)

Construct model's ocean.
Inputs:
    grid    <AbstractGrid> model grid
    u       <Real> ocean x-velocity for each grid line
    v       <Real> ocean y-velocity for each grid line
    temp    <Real> temperature at ocean/ice interface per grid cell
Output: 
        Ocean with constant velocity and temperature on each grid line.
"""
Ocean{FT}(
    grid::AbstractGrid,
    u,
    v,
    temp,
) where {FT <: AbstractFloat} =
    Ocean{FT}(
        fill(FT(u), grid.Nx + 1, grid.Ny + 1),
        fill(FT(v), grid.Nx + 1, grid.Ny + 1),
        fill(FT(temp), grid.Nx + 1, grid.Ny + 1),
    )
