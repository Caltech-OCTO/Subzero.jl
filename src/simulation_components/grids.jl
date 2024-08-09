export AbstractGrid, RegRectilinearGrid

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
    YourGridType{FT} <: AbstractGrid{FT}

Each simulation run with Subzero.jl must be run on a grid. A grid defines the points at
which the ocean and atmosphere hold velocity data, as well as various other tracers at grid
points. The user must choose an existing subtype of `AbstractGrid` or implement a new subtype to
create a simulation. 

Each concrete implementation of `AbstractGrid` must have a field called `floe_locations`
that is a matrix of `CellFloes`[@ref] for each grid cell...

## _API_
The following methods must be implemented for all subtypes:

"""
abstract type AbstractGrid{FT<:AbstractFloat} end

# Concrete type of AbstractGrid - see documentation below
struct RegRectilinearGrid{FT}<:AbstractGrid{FT}
    Nx::Int
    Ny::Int
    Δx::FT
    Δy::FT
    x0::FT
    xf::FT
    y0::FT
    yf::FT
    floe_locations::Matrix{CellFloes{FT}}
end

"""
    RegRectilinearGrid{FT} <: AbstractGrid{FT}

A concrete implementation of an AbstractGrid that represents a tessellation of 2-dimensional
Euclidean space into `n`-by-`m` congruent rectangles. Fields that hold float data are of
type `FT`, a concrete implementation of `AbstractFloat`.

## _Fields_
- `Nx::Int`: number of grid cells in the x-direction
- `Ny::Int`: number of grid cells in the y-direction
- `x0::FT`: value of first x grid line
- `xf::FT`: value of final x grid line
- `y0::FT`: value of first y grid line
- `yf::FT`: value of final y grid line
- `Δx::FT`: grid cell width
- `Δy::FT`: grid cell height
- `floe_locations::Matrix{CellFloes}`: `n + 1` by `m + 1` matrix of [`CellFloes`](@ref)

If Δx doesn't evenly divide x length (lu-lx) or Δy doesn't evenly divide y
length (uy-ly) you won't get full size grid. The grid will be "trimmed" to
the nearest full grid square in both directions.

```jldoctest
julia> grid = RegRectilinearGrid(; x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5, Nx = 20, Ny = 10)
RegRectilinearGrid{Float64}
⊢x extent (0.0 to 500000.0) with 20 grid cells of size 25000.0 m
∟y extent (0.0 to 500000.0) with 10 grid cells of size 50000.0 m
```

```jldoctest
julia> grid = RegRectilinearGrid(Float32; x0 = -5e5, xf = 5e5, y0 = 0.0, yf = 5e5, Δx = 5e4, Δy = 5e4)
RegRectilinearGrid{Float32}
⊢x extent (-500000.0 to 500000.0) with 20 grid cells of size 50000.0 m
∟y extent (0.0 to 500000.0) with 10 grid cells of size 50000.0 m
```

```jldoctest
julia> grid = RegRectilinearGrid(; x0 = -5e5, xf = 5e5, y0 = 0.0, yf = 5e5, Δx = 5e4, Ny = 10)
ERROR: ArgumentError: To create a RegRectilinearGrid, either provide Δx and Δy OR Nx and Ny.
[...]
```

"""
function RegRectilinearGrid(::Type{FT} = Float64;
    x0, xf, y0, yf,
    Nx = nothing, Ny = nothing,
    Δx = nothing, Δy = nothing,
) where FT <: AbstractFloat
    # ensure x and y endpoints are properly ordered
    if x0 > xf
        @warn "$x0 can't be larger than $xf. Swapping the order of the grid x-endpoints"
        xf, x0 = x0, xf
    end
    if y0 > yf
        @warn "$y0 can't be larger than $yf. Swapping the order of the grid y-endpoints"
        yf, y0 = y0, yf
    end
    # determine number and size of gridpoints based on provided arguments
    if !isnothing(Nx) && !isnothing(Ny)  # determine size of grid cells based on number
        Δx = (xf - x0) / Nx
        Δy = (yf - y0) / Ny
    elseif !isnothing(Δx) && !isnothing(Δy)  # determine number of grid cells based on size
        Nx = floor(Int, (xf - x0) / Δx)
        Ny = floor(Int, (yf - y0) / Δy)
        # note: round grid size if requested cell size doesn't go evenly into endpoints
        xf = x0 + Nx * Δx
        yf = y0 + Ny * Δy
    else
        throw(ArgumentError("To create a RegRectilinearGrid, either provide Δx and Δy OR Nx and Ny."))
    end
    # create cell floe list
    floe_locations = [CellFloes{FT}() for _ in 1:(Nx + 1), _ in 1:(Ny + 1)]
    # create RegRectilinearGrid object
    return RegRectilinearGrid{FT}(Nx, Ny, Δx, Δy, x0, xf, y0, yf, floe_locations)
end

# Pretty printing for RegRectilinearGrid showing key dimensions
function Base.show(io::IO, grid::RegRectilinearGrid{FT}) where FT
    overall_summary = "RegRectilinearGrid{$FT}"
    x_summary = "x extent ($(grid.x0) to $(grid.xf)) with $(grid.Nx) grid cells of size $(grid.Δx) m"
    y_summary = "y extent ($(grid.y0) to $(grid.yf)) with $(grid.Ny) grid cells of size $(grid.Δy) m"
    print(io, overall_summary, "\n",
        "⊢", x_summary, "\n",
        "∟", y_summary)
end
