"""
Structs and functions used to define a Subzero grids
"""

"""
    CellFloes{FT<:AbstractFloat}

Struct that tracks which floes are within given cell, as well as the translation
vector needed to move floe each from current postion into cell if it is in cell
due to periodic boundaries. Each index in floeidx is the index of a floe within
the cell and the Δx and Δy with the same index are that floe's translation
vector.
Note: the floeidx is the index of grid cells centered on grid lines, not on the
grid cells defined by the regular, rectilinear grid. 
"""
struct CellFloes{FT<:AbstractFloat}
    floeidx::Vector{Int}
    Δx::Vector{FT}
    Δy::Vector{FT}
end

"""
    CellFloes{FT}()

Constructs an CellFloes object with empty lists for fields.
"""
CellFloes{FT}() where {FT} = CellFloes{FT}(
    Vector{Int}(),
    Vector{FT}(),
    Vector{FT}()
)
"""
    empty!(cell::CellFloes)

Empties the vectors within a CellFloes object
"""
function Base.empty!(cell::CellFloes)
    empty!(cell.floeidx)
    empty!(cell.Δx)
    empty!(cell.Δy)
    return
end

"""
    AbstractGrid

An abstract type for the grid that model will be simulated on.
Affects calculation on the grid.
"""
abstract type AbstractGrid{FT<:AbstractFloat} end

"""
    RegRectilinearGrid{FT<:AbstractFloat}<:AbstractGrid

Tessellation of 2-dimensional Euclidean space into n-by-m congruent rectangles.
- Nx: number of grid cells in the x-direction
- Ny: number of grid cells in the y-direction
- x0: value of first x grid line
- xf: value of final x grid line
- y0: value of first y grid line
- yf: value of final y grid line
- Δx: grid cell width
- Δy: grid cell height
- floe_locations: matrix of CellFloes, which keep track of which floes are in
    each cell
"""
struct RegRectilinearGrid{FT}<:AbstractGrid{FT}
    Nx::Int
    Ny::Int
    x0::FT
    xf::FT
    y0::FT
    yf::FT
    Δx::FT
    Δy::FT
    floe_locations::Matrix{CellFloes{FT}}

    function RegRectilinearGrid{FT}(
        Nx,
        Ny,
        x0,
        xf,
        y0,
        yf,
        Δx,
        Δy,
        floe_locations,
    ) where {FT}
        if size(floe_locations) != (Nx + 1, Ny + 1)
            throw(ArgumentError("Floe location matrix needs to be the same \
                dimensions as grid lines."))
        end
        if (xf - x0) / Δx != Nx
            throw(ArgumentError("X-grid extent and grid cell width don't match \
                with number of grid cells in x-direction."))
        end
        if (yf - y0) / Δy != Ny
            throw(ArgumentError("Y-grid extent and grid cell height don't match \
                with number of grid cells in y-direction."))
        end
        return new{FT}(Nx, Ny, x0, xf, y0, yf, Δx, Δy, floe_locations)
    end
end

"""
    RegRectilinearGrid(::Type{FT}, args...)

A float type FT can be provided as the first argument of any RegRectilinearGrid
constructor. A RegRectilinearGrid of type FT will be created by passing all
other arguments to the correct constructor. 
"""
RegRectilinearGrid(::Type{FT}, args...) where {FT <: AbstractFloat} =
    RegRectilinearGrid{FT}(args...)

"""
    RegRectilinearGrid(args...)

If a type isn't specified, RegRectilinearGrid will be of type Float64 and the
correct constructor will be called with all other arguments.
"""
RegRectilinearGrid(args...) = RegRectilinearGrid{Float64}(args...)

"""
    RegRectilinearGrid{FT}(
        xbounds::Tuple,
        ybounds::Tuple,
        Δx,
        Δy,
    )

Construct a RegRectilinearGrid for model given bounds for grid x and y and grid
cell dimensions in meters.
Inputs:
    xbounds  <Tuple{Real, Real}> bound of grid x-direction in form (left, right)
    ybounds  <Tuple{Real, Real}> bound of grid y-direction in form (bottom, top)
    Δx       <Real> length/height of grid cells in x-direction
    Δy       <Real> length/height of grid cells in y-direction
Output: 
    RegRectilinearGrid from lx to ux and height from ly to uy with grid squares
    of size Δx by Δy
Warning:
    If Δx doesn't evenly divide x length (lu-lx) or Δy doesn't evenly divide y
    length (uy-ly) you won't get full size grid. The grid will be "trimmed" to
    the nearest full grid square in both directions.
"""
function RegRectilinearGrid{FT}(
    xbounds::Tuple,
    ybounds::Tuple,
    Δx,
    Δy,
) where {FT <: AbstractFloat}
    Nx = floor(Int, (xbounds[2] - xbounds[1]) / Δx)
    Ny = floor(Int, (ybounds[2] - ybounds[1]) / Δy)
    return RegRectilinearGrid{FT}(
        Nx,
        Ny,
        xbounds[1],
        xbounds[1] + Nx * Δx,
        ybounds[1],
        ybounds[1] + Ny * Δy,
        Δx,
        Δy,
        [CellFloes{FT}() for i in 1:(Nx + 1), j in 1:(Ny + 1)],
    )
end

get_grid_lines(g0, gf, Δg) = g0:Δg:gf
get_grid_centers(g0, gf, Δg) = (g0 + Δg/2):Δg:(gf - Δg/2) 

"""
    RegRectilinearGrid{FT}(
        Nx,
        Ny,
        xbounds::Tuple{Real, Real},
        ybounds::Tuple{Real, Real},
    ) where {FT <: AbstractFloat}

Construct a RegRectilinearGrid for model given bounds for grid x and y and the
number of grid cells in both the x and y direction.
Inputs:
    Nx       <Int> number of grid cells in the x-direction
    Ny       <Int> number of grid cells in the y-direction
    xbounds  <Tuple{Real, Real}> bound of grid x-direction in form (left, right)
    ybounds  <Tuple{Real, Real}> bound of grid y-direction in form (bottom, top)
Output: 
    RegRectilinearGrid with width and height determined by xbound and ybounds
    and the number of grid cells in the x-direction and y-direction determined
    by dims.
"""
RegRectilinearGrid{FT}(
    Nx::Int,
    Ny::Int,
    xbounds::Tuple,
    ybounds::Tuple,
) where {FT <: AbstractFloat} = 
    RegRectilinearGrid{FT}(
        Nx,
        Ny,
        xbounds[1],
        xbounds[2],
        ybounds[1],
        ybounds[2],
        (xbounds[2] - xbounds[1]) / Nx,
        (ybounds[2] - ybounds[1]) / Ny,
        [CellFloes{FT}() for i in 1:(Nx + 1), j in 1:(Ny + 1)],
    )
