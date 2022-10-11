"""
Structs and functions to create and run a Subzero simulation
"""

"""
    CoarseGridData{FT<:AbstractFloat}

Data averaged over a coarse grid covering the domain. Each field is the size of the coarse grid, with one value per coarse grid square, which is the average of mutliple grid squares in finer resolution grid that the simulation is run using. This data is used for plotting and is saved to be able to examine model progress and data. If all data is desired on the same resolution as the model, the coarse grid should be equal to the model grid and therefore the size of each of these fields will be the same as the model grid's dimensions. Each field must have the same dimensions to be a valid struct.
"""
struct CoarseGridData{FT<:AbstractFloat}
    u::Matrix{FT}
    v::Matrix{FT}
    du::Matrix{FT}
    dv::Matrix{FT}
    stress::Matrix{FT}
    stressxx::Matrix{FT}
    stressxy::Matrix{FT}
    stressyx::Matrix{FT}
    stressyy::Matrix{FT}
    strainux::Matrix{FT}
    strainuy::Matrix{FT}
    strainvx::Matrix{FT}
    strainvy::Matrix{FT}
    frac_ice::Matrix{FT}
    over::Matrix{FT}
    mtot::Matrix{FT}
    area::Matrix{FT}
    height::Matrix{FT}

    CoarseGridData(u, v, du, dv, stress, stressxx, stressxy, stressyx, stressyy,
                 strainux, strainuy, strainvx, strainvy, frac_ice, over, mtot, area, height) = 
                 (size(u) == size(v) == size(du) == size(dv) == size(stress) == size(stressxx) == size(stressxy) == size(stressyx) == 
                 size(stressyy) == size(strainux) ==  size(strainuy) ==
                 size(strainvx) == size(strainvy) == size(frac_ice) ==
                 size(over) == size(mtot) == size(area) == size(height)) ?
                 new{eltype(u)}(u, v, du, dv, stress, stressxx, stressxy, stressyx, stressyy, strainux, strainuy, strainvx, strainvy, frac_ice, over, mtot, area, height) :
                 throw(ArgumentError("Size of coarse grid data matrices are not uniform between all of the fields."))
end

"""
    CoarseGridData(coarse_nx, coarse_ny, t::Type{T} = Float64)

Initialize all fields to correct size matricies where all values are 0s.

Inputs:
        coarse_nx   <Int> number of grid squares in the x-direction within the 
                        coarse grid
        coarse_ny   <Int> number of grid squares in the y-direction within the 
                        coarse grid
        t           <Type> datatype to convert grid fields - must be a Float!
Outputs:
        CoarseGridData object where each field is a matrix of dimensions coarse_nx by coarse_ny and all elements of each matrix are 0 of type T. 
"""
function CoarseGridData(coarse_nx, coarse_ny, t::Type{T} = Float64) where T
    u = zeros(T, coarse_nx, coarse_ny)
    v = zeros(T, coarse_nx, coarse_ny)
    du = zeros(T, coarse_nx, coarse_ny)
    dv = zeros(T, coarse_nx, coarse_ny)
    stress = zeros(T, coarse_nx, coarse_ny)
    stressxx = zeros(T, coarse_nx, coarse_ny)
    stressxy = zeros(T, coarse_nx, coarse_ny)
    stressyx = zeros(T, coarse_nx, coarse_ny)
    stressyy = zeros(T, coarse_nx, coarse_ny)
    strainux = zeros(T, coarse_nx, coarse_ny)
    strainuy = zeros(T, coarse_nx, coarse_ny)
    strainvx = zeros(T, coarse_nx, coarse_ny)
    strainvy = zeros(T, coarse_nx, coarse_ny)
    frac_ice = zeros(T, coarse_nx, coarse_ny)
    over = zeros(T, coarse_nx, coarse_ny)
    mtot = zeros(T, coarse_nx, coarse_ny)
    area = zeros(T, coarse_nx, coarse_ny)
    height = zeros(T, coarse_nx, coarse_ny)
    return CoarseGridData(u, v, du, dv, stress, stressxx, stressxy, stressyx,
                          stressyy, strainux, strainuy, strainvx, strainvy, frac_ice, over, mtot, area, height)
end

"""
    CoarseGridData(grid::Grid, t::Type{T} = Float64)

Initialize all fields to correct size matricies where all values are 0s.

Inputs:
        coarse_nx   <Int> number of grid squares in the x-direction within the 
                        coarse grid
        coarse_ny   <Int> number of grid squares in the y-direction within the 
                        coarse grid
        t           <Type> datatype to convert grid fields - must be a Float!
Outputs:
        CoarseGridData object where each field is a matrix of dimensions coarse_nx by coarse_ny and all elements of each matrix are 0 of type T. 
"""
function CoarseGridData(grid::Grid, t::Type{T} = Float64) where T
    return CoarseGridData(grid.dims[2], grid.dims[1], T)
end

"""
    Simulation{FT<:AbstractFloat, DT<:AbstractDomain{FT}}

Simulation which holds a model and parameters needed for running the simulation. Simulation requires a model, a coarse grid, a coarse grid data struct, and a figure. The figure can be initialized using setup_plot. The rest of the simulation values are optional. These fields and their default values are as follows:  the size of a timestep in seconds Δt (10), the total number of timesteps in the simulation nΔt (7500), the output frequency of floe and data on the coarse grid in timesteps nΔtout (150),  timesteps between saving images Δtpics (150),  timesteps between floe simplicaiton  Δtsimp (20), timesteps betwen thermodynamic floe creation Δtpack (500), timesteps between updating ocean forcing  Δtocn (10). There are also flags that control simulation behavior. These flags are AVERAGE (average coarse grid data in time), COLLISION (enable floe collisions), CORNERS (floe corners can break), FRACTURES (floes can fracture), KEEPMIN (small floes don't dissolve), PACKING (floe packing enabled), RAFTING (floe rafting enabled), RIDGING (floe ridging enabled), and WELDING (floe welding enabled). All are false by default.
"""
@kwdef struct Simulation{FT<:AbstractFloat, DT<:AbstractDomain{FT}}
    # Objects ------------------------------------------------------------------
    model::Model{FT, DT}            # Model to simulate
    cgrid::Grid{FT}                 # Coarse grid used to average and save 
                                    # output values
    cgrid_data::CoarseGridData{FT}  # Saved values from coarse grid calculations
    fig::Plots.Plot                 # Figure used for basic plotting of floes - 
                                    # use setup_plot function to generate!
    # Timesteps ----------------------------------------------------------------
    Δt::Int = 10                    # Simulation timestep (seconds)
    nΔt::Int = 7500                 # Total timesteps simulation runs for
    nΔtout::Int = 150               # Output frequency in timesteps to save
                                    # cgrid_data and floe data
    Δtpics::Int = 150               # Timesteps between saving images
    Δtsimp::Int = 20                # Timesteps between floe simplification
    Δtpack::Int = 500               # Timesteps between thermodynamic floe 
                                    # creation
    Δtocn::Int = 10                 # Timesteps between updating ocean forces
    # Flags --------------------------------------------------------------------
    AVERAGE::Bool = false           # If true, average coarse grid data in time
    COLLISION::Bool = false         # If true, collisions are enabled for floes
    CORNERS::Bool = false           # If true, corners of floes can break
    FRACTURES::Bool = false         # If true, fracturing of floes is enabled
    KEEPMIN::Bool = false           # If true, retain small floes that would 
                                    # normally dissolve
    PACKING::Bool = false           # If true, floe packing is enabled
    RAFTING::Bool = false           # If true, floe rafting is enabled
    RIDGING::Bool = false           # If true, floe ridging is enabled
    WELDING::Bool = false           # If true, floe welding is enabled
end

# function calc_ocean_coupling(grid, floe_arr, topography_arr, coarse_nx, coarse_ny, PERIODIC)
#     live_floes = filter(floe->floe.alive, floe_arr)
#     topography_poly = LG.MultiPolygon([translate(poly.coords, centroid) for
#                                        poly in topography_arr])
#     # TODO: Periodic section!

#     # create coarse grid
# end

"""
    domain_coords(domain::RectangleDomain)
Inputs:
        domain<RectangleDomain>
Output:
        RingVec coordinates for edges of rectangular domain based off of boundary values
"""
function domain_coords(domain::RectangleDomain)
    northval = domain.north.val
    southval = domain.south.val
    eastval = domain.east.val
    westval = domain.west.val
    coords = [[westval, northval], [westval, southval],
              [eastval, southval], [eastval, northval],
              [westval, northval]]
    return coords
end #TODO: Might not need!

function run!(simulation)
    plot_sim(simulation.model, simulation.fig, 1)
end