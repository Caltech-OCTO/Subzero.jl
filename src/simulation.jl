"""
Structs and functions to create and run a Subzero simulation
"""


#= function calc_eulerian(grid, floe_arr, topography_arr, coarse_nx, coarse_ny, PERIODIC)
    live_floes = filter(floe->floe.alive, floe_arr)
    topography_poly = LG.MultiPolygon([translate(poly.coords, centroid) for
                                       poly in topography_arr])
    # TODO: Periodic section!

    # create coarse grid
end =#

struct EulerianGrid{FT<:AbstractFloat}
    coarse_nx ::Int
    coarse_ny ::Int
    xlines::Matrix{FT}
    ylines::Matrix{FT}
    xgrid::Matrix{FT}
    ygrid::Matrix{FT}

    EulerianGrid(xgrid, ygrid) =  (size(xgrid) == size(ygrid)) ?
        new{eltype(xgrid)}(xgrid, ygrid) :
        throw(ArgumentError("Grid dimensions don't match within EularianGrid."))
end

#= function EulerianGrid(domain::CircleDomain, coarse_nx::Int, coarse_ny::Int)
    nval = domain.centroid[2] + domain.radius
    sval = domain.centroid[2] - domain.radius
    eval = domain.centroid[1] + domain.radius
    wval = domain.centroid[1] - domain.radius
    e_xlines = wval:(eval - wval)/coarse_nx:eval
    e_ylines = sval:(nval - sval)/coarse_ny:nval
    # do meshgrid
end
 =#
#= function EulerianGrid(domain::RectangleDomain, coarse_nx::Int, coarse_ny::Int)
    nval = domain.north.val
    sval = domain.south.val
    eval = domain.east.val
    wval = domain.west.val
    e_xlines = wval:(eval - wval)/coarse_nx:eval
    e_ylines = sval:(nval - sval)/coarse_ny:nval
    # do meshgrid
end =#

struct EulerianData{FT<:AbstractFloat}
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

    EularianData(u, v, du, dv, stress, stressxx, stressxy, stressyx, stressyy,
                 strainux, strainuy, strainvx, strainvy, frac_ice, over, mtot, area, height) = 
                 (size(u) == size(v) == size(du) == size(dv) == size(stress) == size(stressxx) == size(stressxy) == size(stressyx) == 
                 ize(stressyy) == size(strainux) ==  size(strainuy) ==
                 size(strainvx) == size(strainvy) == size(frac_ice) ==
                 size(over) == size(mtot) == size(area) == size(height)) ?
                 new{eltype(u)}(u, v, du, dv, stress, stressxx, stressxy, stressyx, stressyy, strainux, strainuy, strainvx, strainvy, frac_ice, over, mtot, area, height) :
                 throw(ArgumentError("Size of Eularian data matrices are not uniform between all fields."))
end

function EulerianData(coarse_nx, coarse_ny, t::Type{T} = Float64) where T
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
end

struct Simulation
    model::Model
    PERIODIC::Bool
    RIDGING::Bool
    FRACTURES::Bool
    PACKING::Bool
    WELDING::Bool
    CORNERS::Bool
    COLLISION::Bool
    RAFTING::Bool
    AVERAGE::Bool
    KEEP_MIN::Bool
end