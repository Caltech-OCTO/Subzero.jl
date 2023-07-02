"""
Structs and functions used to define a Subzero atmosphere
"""


"""
Atmos velocities in the x-direction (u) and y-direction (v). u and v should
match the size of the corresponding model grid so that there is one x and y
velocity value for each grid cell. Atmos also needs temperature at the
atmosphere/ice interface in each grid cell. Model cannot be constructed if size
of atmos fields and grid do not match.
"""
struct Atmos{FT<:AbstractFloat}
    u::Matrix{FT}
    v::Matrix{FT}
    temp::Matrix{FT}
    function Atmos{FT}(u, v, temp) where {FT <: AbstractFloat}
        if !(size(u) == size(v))
            throw(ArgumentError("One or more of the atmosphere vector fields \
             aren't the same dimension."))
        end
        new{FT}(u, v, temp)
    end
end

"""
    Atmos(::Type{FT}, args...)

A float type FT can be provided as the first argument of any Atmos constructor.
An Atmos of type FT will be created by passing all other arguments to the
correct constructor. 
"""
Atmos(::Type{FT}, args...) where {FT <: AbstractFloat} =
    Atmos{FT}(args...)

"""
    Atmos(args...)

If a type isn't specified, Atmos will be of type Float64 and the correct
constructor will be called with all other arguments.
"""
Atmos(args...) = Atmos{Float64}(args...)

"""
    Atmos{FT}(grid, u, v)

Construct model atmosphere of type FT.
Inputs:
    grid        <AbstractGrid> model's grid 
    u           <Real> Atmos x-velocity for each grid cell
    v           <Real> Atmos y-velocity for each grid cell
    temp        <Real> temperature at atmopshere/ice interface per grid cell
Output: 
    Atmosphere of type FT with constant velocity and temperature over domain.
"""
Atmos{FT}(
    grid::AbstractGrid,
    u,
    v,
    temp,
) where {FT <: AbstractFloat} = 
    return Atmos{FT}(
        fill(FT(u), grid.Nx + 1, grid.Ny + 1),
        fill(FT(v),  grid.Nx + 1, grid.Ny + 1),
        fill(FT(temp),  grid.Nx + 1, grid.Ny + 1),
    )
