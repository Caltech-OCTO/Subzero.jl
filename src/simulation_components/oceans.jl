export Ocean

# See CellFloes documentation below
struct CellStresses{FT<:AbstractFloat}
    τx::Vector{FT}
    τy::Vector{FT}
    npoints::Vector{Int}
end

"""
    CellStresses{FT}

Struct to collect stress from ice floes on ocean grid cells. One `CellStresses`
corresponds to one grid cell and is "linked" to a corresponding [`CellFloes`](@ref) object.
The `CellFloes` object records which floes are in a given grid cell and the `CellStresses`
aggregates the total stress felt on the ocean from the ice (on each of the floe's sub-floe
points generated by a subtype of [`AbstractSubFloePointsGenerator`](@ref)) in a given grid
cell. This is needed ONLY for two-way coupling. Due to the prevalence of periodic
boundaries, the grid cell represented by a `CellStresses` object are centered on grid
points (translated by `Δx/2` and `Δy/2` in the x and y directions for a
[`RegRectilinearGrid`](@ref)), rather than on the grid cells defined by the grid object
itself. A `CellStresses` object holds a list of stresses on the ocean where each element is
from a single floe within grid cell. The ith element of the `τx` and `τy` fields are measures
of the total stress from the ith floe in the cell and the `npoints` field is a measure of
how many of the ith floe's sub-floe points contributed to the running total of the stress.
Again, each element in the list corresponds to one floe, which is denoted in the
corresponding [`CellFloes`](@ref) matrix within the grid. 

##  _Fields_
- `τx::Vector{FT}`: list of total x-stress caused by each floe in a cell on the ocean
- `τy::Vector{FT}`: list of total y-stress caused by each floe in a cell on the ocean
- `npoints::Vector{Int}`: list of total number of sub-floe points contributing to each individial floe's stress

Here is how to construct a `CellStresses` object:

    CellStresses([FT = Float64]; τx = nothing, τy = nothing, npoints = nothing)

## _Positional arguments_
- $FT_DEF
## _Keyword arguments_
- `τx::Vector{FT}`: list of total x-stress caused by each floe in a cell on the ocean
- `τy::Vector{FT}`: list of total y-stress caused by each floe in a cell on the ocean
- `npoints::Vector{Int}`: list of total number of sub-floe points contributing to each individial floe's stress

!!! note
    If no keyword arguments are provide by the user, an `CellStresses` object with empty
    fields will be created. This is the **standard useage** of these objects and they are added to
    during the coupling step. If keyword arguments are provided, then all three must be provided
    and each vector must be the same size.

"""
function CellStresses(::Type{FT} = Float64; τx = nothing, τy = nothing, npoints = nothing) where FT
    if isnothing(τx) || isnothing(τy) || isnothing(npoints)
        τx = Vector{FT}()
        τy = Vector{FT}()
        npoints = Vector{Int}()
    else
        @assert length(τx) == length(τy) == length(npoints) "A CellStresses object requires that all fields be vectors of the same length."
    end
    return CellStresses{FT}(τx, τy, npoints)
end

# Empties each of the three vectors (`τx`, `τy`, and `npoints`) within an CellStresses object
function Base.empty!(cell::CellStresses)
    empty!(cell.τx)
    empty!(cell.τy)
    empty!(cell.npoints)
    return
end

const OCEAN_TEMP_STR = "Ocean temperatures are above the range for freezing. The thermodynamics aren't currently setup for these conditions."

# See documentation below
struct Ocean{FT<:AbstractFloat}
    u::Matrix{FT}
    v::Matrix{FT}
    temp::Matrix{FT}
    hflx_factor::Matrix{FT}
    scells::Matrix{CellStresses{FT}}
    τx::Matrix{FT}
    τy::Matrix{FT}
    si_frac::Matrix{FT}
    dissolved::Matrix{FT}

    function Ocean{FT}(u, v, temp, hflx_factor, scells, τx, τy, si_frac, dissolved) where FT
        if !all(-3 .<= temp .<= 0)
            @warn OCEAN_TEMP_STR
        end
        if !(size(u) == size(v) == size(τx) == size(τy))
            throw(ArgumentError("One or more of the ocean vector fields aren't \
                the same dimension."))
        end
        if !(size(temp) == size(hflx_factor) == size(si_frac) == size(dissolved))
            throw(ArgumentError("One or more of the ocean tracer fields aren't \
                the same dimension."))
        end
        new{FT}(u, v, temp, hflx_factor, scells, τx, τy, si_frac, dissolved)
    end
end

const U_DEF = "`u::Union{Real, Matrix{Real}}`: user can either provide a single u-velocity value \
(a constant field the size of the `grid` will be created) or provide a matrix field of  u-velocity values"
const V_DEF = "`v::Union{Real, Matrix{Real}}`: user can either provide a single v-velocity value \
(a constant field the size of the `grid` will be created) or provide a matrix field of v-velocity values"
const TEMP_DEF = "`temp::Union{Real, Matrix{Real}}`: user can either provide a single temperature value \
(a constant field the size of the `grid` will be created) or provide a matrix field of temperature values"

"""
    Ocean{FT}

Each simulation needs ocean values on the grid to drive the ice floes. The `Ocean` struct
holds the needed fields to perform coupling, including vector fields (like u-velocity,
v-velocity, x-stress, and y-stress) and tracer fields (like temperature, heatflux, and 
dissolved ice mass). We also hold and calculate the sea ice fraction in any given grid cell.

Right now, all values are stored on grid points as defined by an [`AbstractRectilinearGrid`](@ref).
This should eventually change to a c-grid. If interested in this change, see issue [#30](@ref).
This will require a big change in behavior, but is an important change for consistency and
to match with Oceananigans for two-way coupling.

Coupling beyond velocity (for example thermodynamically) is not very well developed. If you
are interested in working on this see issue [#9](@ref). The `hflx_factor` field is related
to thermodynamic coupling.

## _Fields_
- `u::Matrix{FT}`: ocean u-velocities in the x-direction on each grid point
- `v::Matrix{FT}`: ocean v-velocities in the y-direction on each grid point
- `temp::Matrix{FT}`: ocean temperature on each grid point
- `hflx_factor::Matrix{FT}`: factor to calculate the ocean-atmosphere heat flux for a grid cell
- `scells::Matrix{CellStresses}`: used to accumulate stresses on the ocean per grid cell from floes in cell (see [`CellStresses`](@ref))
- `τx::Matrix{FT}`: stress on the ocean in the x-direction on each grid point
- `τy::Matrix{FT}`: stress on the ocean in the y-direction on each grid point
- `si_frac::Matrix{FT}`: fraction of area in each grid cell (centered on grid points) that is covered in sea-ice (between 0-1)
- `dissolved::Matrix{FT}`: total mass from ice floes dissolved in each grid cell (centered on grid points)

!!! note
    If a periodic boundary is used in the domain, the last grid cell in that direction will
    not be used as it is equivalent to the first grid cell. Thus, for all of these fields, the
    first and last value in the x and/or y direction should be equal if the east-west or
    north-south boundary pair are periodic respectively. In the future, it might be worth having
    the ocean re-sizing to not included this repeated point dispatching off of the boundary types.

!!! note
    For each of the fields, the rows represent the `x`-values of the grid, while the columns
    represent the `y`-values of the grid. This makes it easy to index into the grid for a
    point `(x, y)`, although it does mean that the fields provided don't "look" like the grid
    if plotted directly and instead need to be transposed. 

Here is how to construct an `Ocean`:

    Ocean([FT = Float64]; u, v, temp, grid = nothing)

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- $U_DEF
- $V_DEF
- $TEMP_DEF
- $GRID_DEF

!!! note
    The user MUST provide `u`, `v`, and `temp`, although they have the option of these
    being `Real` values or `Matrix{<:Real}`. If the user chooses to provide `Real` values for any
    of those fields, then a `grid` IS required to create a constant field of the correct size.
    If user chooses to provide a field, they should be of size `(grid.Nx + 1, grid.Ny + 1)` and
    the user doesn't require a `grid` input.

## _Examples_
- Creating `Ocean` with constant `u`-velocity, `v-velocity`, and `temp`
```jldoctest ocean
julia> grid = RegRectilinearGrid(; x0 = 0.0, xf = 5e5, y0 = 0.0, yf = 5e5, Nx = 20, Ny = 10);

julia> Ocean(; u = 0.5, v = 0.25, temp = 0.0, grid)
Ocean{Float64}
  ⊢Vector fields of dimension (21, 11)
  ⊢Tracer fields of dimension (21, 11)
  ⊢Average u-velocity of: 0.5 m/s
  ⊢Average v-velocity of: 0.25 m/s
  ∟Average temperature of: 0.0 C
```
- Creating `Ocean` with user-defined `u` and `v` and a constant `temp` with `Float32` data
```jldoctest ocean
julia> import Random.Xoshiro;

julia> seeded_rng = Xoshiro(100);

julia> u_field = rand(seeded_rng, grid.Nx + 1, grid.Ny + 1);

julia> v_field = rand(seeded_rng, grid.Nx + 1, grid.Ny + 1);

julia> Ocean(Float32; u = u_field, v = v_field, temp = 0.0, grid)
Ocean{Float32}
  ⊢Vector fields of dimension (21, 11)
  ⊢Tracer fields of dimension (21, 11)
  ⊢Average u-velocity of: 0.4856711 m/s
  ⊢Average v-velocity of: 0.5319979 m/s
  ∟Average temperature of: 0.0 C
```
- Trying to create an `Ocean` with a constant field and NO grid
```jldoctest
julia> Ocean(; u = 0.2, v = 0.1, temp = 0.0)
ERROR: ArgumentError: To create a matrix from the constant value provided, you must provide a grid.
[...]
```
"""
function Ocean(::Type{FT} = Float64; u, v, temp, grid = nothing) where {FT <: AbstractFloat}
    # Create field from provided values
    u_field = _get_val_field(FT, u, grid)
    v_field = _get_val_field(FT, v, grid)
    temp_field = _get_val_field(FT, temp, grid)
    # Create rest of fields using the size of the temperature field
    Nx, Ny = size(temp_field)
    hflx_factor = zeros(FT, Nx, Ny) # heat flux
    scells = [CellStresses(FT) for _ in 1:Nx, _ in 1:Ny]
    τx = zeros(FT, Nx, Ny)
    τy = zeros(FT, Nx, Ny)
    si_frac = zeros(FT, Nx, Ny)
    dissolved = zeros(FT, Nx, Ny)
    # Create Ocean from fields
    return Ocean{FT}(u_field, v_field, temp_field, hflx_factor, scells, τx, τy, si_frac, dissolved)
end

# If matrix of values and NO grid is provided, convert to type FT and return
_get_val_field(::Type{FT}, v::AbstractArray, grid) where FT = FT.(v)
# If single value and grid is provided, make constant matrix filled with provided value converted to type FT anf return
_get_val_field(::Type{FT}, v::Real, grid::AbstractRectilinearGrid) where FT = fill(FT(v), grid.Nx + 1, grid.Ny + 1)
# If single value and NO grid is provided, throw an error
_get_val_field(::Type, ::Real, ::Nothing) = throw(ArgumentError("To create a matrix from the constant value provided, you must provide a grid."))
_get_val_field(t, v, g) = throw(ArgumentError("Incorrect inputs to Ocean or Atmos constructor."))

# Pretty printing for Ocean showing key dimensions
function Base.show(io::IO, ocean::Ocean{FT}) where FT
    overall_summary = "Ocean{$FT}"
    vector_summary = "Vector fields of dimension $(size(ocean.u))"
    avg_ocean_u = "Average u-velocity of: $(mean(ocean.u)) m/s"
    avg_ocean_v = "Average v-velocity of: $(mean(ocean.v)) m/s"
    tracer_summary = "Tracer fields of dimension $(size(ocean.temp))"
    avg_ocean_temp = "Average temperature of: $(mean(ocean.temp)) C"
    print(io, overall_summary, "\n",
        "  ⊢", vector_summary, "\n",
        "  ⊢", tracer_summary, "\n",
        "  ⊢", avg_ocean_u, "\n",
        "  ⊢", avg_ocean_v, "\n",
        "  ∟", avg_ocean_temp)
end