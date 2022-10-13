"""
Structs and functions to calculate and write output from the simulation
"""

abstract type AbstractOutputWriter{FT <: AbstractFloat} end

struct GridOutputWriter{FT}<:AbstractOutputWriter{FT}
    nΔtout::Int
end

struct FloeOutputWriter{FT}<:AbstractOutputWriter{FT}
    nΔtout::Int
end


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

nΔtout::Int = 150               # Output frequency in timesteps to save
                                # cgrid_data and floe data
Δtpics::Int = 150               # Timesteps between saving images

AVERAGE::Bool = false           # If true, average coarse grid data in time

# have abstract output writer 

# grid output writer and floe output writer

# each output writer has the same time frequency and same grid size

# ways to get/calculate each field that is within the output writer

