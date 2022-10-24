"""
Structs and functions to calculate and write output from the simulation
"""

"""
    AbstractOutputWriter

An abstract type for output writers that provide data from simulation runs.
"""
abstract type AbstractOutputWriter end


"""
    GridOutputWriter{FT<:AbstractFloat, ST<:AbstractString}<:AbstractOutputWriter

Grid subtype of AbstractOutputWriter that holds information for outputting Eularian data on the grid. This output does not need to be the grid defined for the model. This grid can be coarser, or more fine as defined by the xg, yg, xc, and yc fields. Output on this scale will be saved to the file defined by fn every Δtout timesteps. The output values that are recorded is defined by names, units, and comments.
"""
struct GridOutputWriter{FT<:AbstractFloat, ST<:AbstractString}<:AbstractOutputWriter
    Δtout::Int              # Number of timesteps between grid outputs
    fn::ST                  # Filename for output file
    xg::Vector{FT}          # Grid lines in x-direction for ouput calculations
    yg::Vector{FT}          # Grid lines in y-direction for ouput calculations
    xc::Vector{FT}          # Center-lines on grid cells in the x-direction
    yc::Vector{FT}          # Center-lines on grid cells in the y-direction
    names::Vector{ST}
    units::Vector{ST}
    comments::Vector{ST}

    GridOutputWriter(Δtout, fn, xg, yg, xc, yc, names, units, comments) = 
        (length(xg) > 1 && length(yg) > 1) && length(xg) - 1 == length(xc) && length(yg) - 1 == length(yc) &&
        (length(names) == length(units) == length(comments)) && length(fn)>0 ?
        new{eltype(xg), eltype(names)}(Δtout, fn, xg, yg, xc, yc, names, units, comments) :
        throw(ArgumentError("Output grid lines must have at least one grid square."))
end


"""
    GridOutputWriter(Δtout, fn, grid::Grid, dims)

Create GridOutputWriter for given grid, desired output dimensions to re-grid given grid, and desired number of timesteps, saved to given filename.
Inputs:
        Δtout    <Int> number of timesteps between output
        fn       <String> name of file to save grid data to
        grid     <Grid> original grid, which we are re-gridding
        dims     <(Int, Int)> output grid dimensions - rows -> ny, cols -> nx
Output:
        GridOutputWriter that re-grids given grid to given dimensions, and number of timesteps between saving averaged output on this re-gridded grid to filename.
"""
function GridOutputWriter(Δtout, fn, grid::Grid, dims)
    grid_names = ["u", "v", "du", "dv", "si_frac", "overlap", "mass", "area", "height"]
    grid_units = ["m/s", "m/s", "m/s^2", "m/s^2", "unitless", "m^2", "kg", "m^2", "m"]
    grid_comments = ["Average x-velocity of floes within grid cell", "Average  
                y-velocity of floes within grid cell",
                "Average x-acceleration of floes within grid cell", "Average y-acceleration of floes within grid cell",
                "Fraction of grid cell covered by ice floes", "Average overlap area of floes per floe within grid cell",
                "Mass of floes within grid cell", "Area of floes within grid cell", "Average height of floes within grid cell"]
    lx = grid.xg[1]
    ux = grid.xg[end]
    ly = grid.yg[1]
    uy = grid.yg[end]
    Δx = (ux-lx)/dims[2]
    Δy = (uy-ly)/dims[1]
    xg = collect(lx:Δx:ux) 
    yg = collect(ly:Δy:uy)
    xc = collect(xg[1]+Δx/2:Δx:xg[end]-Δx/2)
    yc = collect(yg[1]+Δy/2:Δy:yg[end]-Δy/2)
    return GridOutputWriter(Δtout, fn, xg, yg, xc, yc, grid_names, grid_units, grid_comments)
end

"""
    FloeOutputWriter{FT<:AbstractFloat, ST<:AbstractString}<:AbstractOutputWriter

Floe subtype of AbstractOutputWriter that holds information for outputting floe information from model throughout simulation. The edges of the grid are defined by xg and yg (the grid lines for the model) for plotting and reconstruction purposes. Output on this scale will be saved to the file defined by fn every Δtout timesteps. The output values that are recorded is defined by names, units, and comments.
"""
struct FloeOutputWriter{FT<:AbstractFloat, ST<:AbstractString}<:AbstractOutputWriter
    Δtout::Int              # Number of timesteps between grid outputs
    fn::ST                  # Filename for output file
    xg::Vector{FT}          # Grid lines in x-direction for plotting
    yg::Vector{FT}          # Grid lines in y-direction for plotting
    names::Vector{ST}
    units::Vector{ST}
    comments::Vector{ST}

    FloeOutputWriter(Δtout, fn, xg, yg, names, units, comments) = 
        length(xg) > 1 && length(yg) > 1 && (length(names) == length(units) == length(comments)) &&
        length(fn) > 0 ?
        new{eltype(xg), eltype(names)}(Δtout, fn, xg, yg, names, units, comments) :
        throw(ArgumentError("Output grid lines must have at least one grid square."))
end

"""
    FloeOutputWriter(Δtout, fn, grid::Grid)

Create FloeOutputWriter for given grid and desired number of timesteps, saved to the given file name.
Inputs:
    Δtout    <Int> number of timesteps between output
    fn       <String> name of file to save grid data to
    grid     <Grid> original grid, which we are re-gridding
Output:
    GridOutputWriter that holds desired floe fields to save to given filename and the number of timesteps between saving this output. Grid extents is saved as metadata for plotting.
"""
function FloeOutputWriter(Δtout, fn, grid::Grid)
    floe_names = ["height", "area", "mass", "moment", "rmax", "α", "xcentroid", 
                  "ycentroid", "u", "v", "ξ", "fxOA", "fyOA", "torqueOA", "p_dxdt", "p_dydt", "p_dudt", "p_dvdt", "p_dξdt", "p_dαdt", "overarea", "alive"]
    floe_units = [ "m", "m^2", "kg", "kg m^2", "m", "rad", "location",
                   "location", "m/s", "m/s", "rad/s", "N", "N", "N m", "m/s", "m/s", "m/s^2", "m/s^2", "rad/s^2", "rad/s", "m^2", "unitless"]
    floe_comments = ["Uniform floe height", "Floe area", "Floe mass",
                     "Floe mass moment of intertia", "Maximum radius within floe", "Floe rotation since starting position", "X-point of centorid", "Y-point of centroid",
                     "Floe x-direction velocity", "Floe y-direction velocity", "Floe angular velocity", "X-directional forces on floe from the ocean and atmosphere",
                     "Y-directional forces on floe from the ocean and atmosphere", "Forque on floes from ocean and atmosphere", "Floe x-velocity from previous time step",
                     "Floe y-velocity from previous timestep", "Floe x-acceleration from previous timestep", "Floe y-acceleration from previous timestep",
                     "Floe angular acceleration from the previous timestep", "Floe angular velocity from previous timestep", "Overlap area of floe with other floes",
                     "Flag if floe is still active in simulation"]
    return FloeOutputWriter(Δtout, fn, grid.xg, grid.yg, floe_names, floe_units, floe_comments)
end


"""
    CoarseGridData{FT<:AbstractFloat}

Data averaged over a coarse grid covering the domain. Each field is the size of the coarse grid, with one value per coarse grid square, which is the average of mutliple grid squares in finer resolution grid that the simulation is run using. This data is used for plotting and is saved to be able to examine model progress and data. If all data is desired on the same resolution as the model, the coarse grid should be equal to the model grid and therefore the size of each of these fields will be the same as the model grid's dimensions. Each field must have the same dimensions to be a valid struct.
"""
struct OutputGridData{FT<:AbstractFloat}
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
    si_frac::Matrix{FT}
    overlap::Matrix{FT}
    mass::Matrix{FT}
    area::Matrix{FT}
    height::Matrix{FT}

    OutputGridData(u, v, du, dv, stress, stressxx, stressxy, stressyx, stressyy,
                 strainux, strainuy, strainvx, strainvy, si_frac, overlap, mtot, area, height) = 
                 (size(u) == size(v) == size(du) == size(dv) == size(stress) == size(stressxx) == size(stressxy) == size(stressyx) == 
                 size(stressyy) == size(strainux) ==  size(strainuy) ==
                 size(strainvx) == size(strainvy) == size(si_frac) ==
                 size(overlap) == size(mtot) == size(area) == size(height)) ?
                 new{eltype(u)}(u, v, du, dv, stress, stressxx, stressxy, stressyx, stressyy, strainux, strainuy, strainvx, strainvy, si_frac, overlap, mtot, area, height) :
                 throw(ArgumentError("Size of coarse grid data matrices are not uniform between all of the fields."))
end

"""
    CoarseGridData(coarse_nx, coarse_ny, t::Type{T} = Float64)

Initialize all fields to correct size matricies where all values are 0s.

Inputs:
        dims        <Vector{Float}> dimensions of output data grid
        t           <Type> datatype to convert grid fields - must be a Float!
Outputs:
    OutputGridData where each field is a matrix of dimensions coarse_nx by coarse_ny and all elements of each matrix are 0 of type T. 
"""
function OutputGridData(dims, t::Type{T} = Float64) where T
    u = zeros(T, dims[1], dims[2])
    v = zeros(T, dims[1], dims[2])
    du = zeros(T, dims[1], dims[2])
    dv = zeros(T, dims[1], dims[2])
    stress = zeros(T, dims[1], dims[2])
    stressxx = zeros(T, dims[1], dims[2])
    stressxy = zeros(T, dims[1], dims[2])
    stressyx = zeros(T, dims[1], dims[2])
    stressyy = zeros(T, dims[1], dims[2])
    strainux = zeros(T, dims[1], dims[2])
    strainuy = zeros(T, dims[1], dims[2])
    strainvx = zeros(T, dims[1], dims[2])
    strainvy = zeros(T, dims[1], dims[2])
    si_frac = zeros(T, dims[1], dims[2])
    overlap = zeros(T, dims[1], dims[2])
    mass = zeros(T, dims[1], dims[2])
    area = zeros(T, dims[1], dims[2])
    height = zeros(T, dims[1], dims[2])
    return OutputGridData(u, v, du, dv, stress, stressxx, stressxy, stressyx,
                          stressyy, strainux, strainuy, strainvx, strainvy, si_frac, overlap, mass, area, height)
end


#AVERAGE::Bool = false           # If true, average coarse grid data in time


function reset_output_grid!(data)
    fill!(data.u, 0.0)
    fill!(data.v, 0.0)
    fill!(data.du, 0.0)
    fill!(data.dv, 0.0)
    fill!(data.stress, 0.0)
    fill!(data.stressxx, 0.0)
    fill!(data.stressxy, 0.0)
    fill!(data.stressyx, 0.0)
    fill!(data.stressyy, 0.0)
    fill!(data.strainux, 0.0)
    fill!(data.strainuy, 0.0)
    fill!(data.strainvx, 0.0)
    fill!(data.strainvy, 0.0)
    fill!(data.si_frac, 0.0)
    fill!(data.overlap, 0.0)
    fill!(data.mass, 0.0)
    fill!(data.area, 0.0)
    fill!(data.height, 0.0)
end

"""
    calc_eulerian_data!(floes, topography, writer, istep)

Calculate floe data averaged on grid defined by GridOutputWriter for current timestep (istep).
Inputs:
        floes       <StructArray{Floe}> array of model's floes
        topography  <StructArray{Topography} array of  model's topography
        writer      <GridOutputWriter> 
        istep       <Int> current simulation timestep
Output:
        Floe data averaged on eularian grid provided. 
"""
function calc_eulerian_data!(floes, topography, writer, istep)
    # Calculate needed values
    live_floes = filter(f -> f.alive == 1, floes)
    Δx = writer.xg[2] - writer.xg[1]
    Δy = writer.yg[2] - writer.yg[1]
    cell_rmax = sqrt(Δx^2 + Δy^2)
    #reset_output_grid!(outdata)
    xgrid, ygrid = grids_from_lines(writer.xc, writer.yc)
    dims = size(xgrid)
    floe_centroids = live_floes.centroid
    floe_rmax = live_floes.rmax

    # Identify floes that potentially overlap each grid square
    potential_interactions = zeros(dims[1], dims[2], length(floes))
    for i in eachindex(floes)
        pint = sqrt.((xgrid .- floe_centroids[i][1]).^2 .+ (ygrid .- floe_centroids[i][2]).^2) .- (floe_rmax[i] + cell_rmax)
        pint[pint .> 0] .= 0
        pint[pint .< 0] .= 1
        potential_interactions[:,:,i] = pint
    end

    ds = NCDataset(joinpath(pwd(), "output", "grid", writer.fn), "a")
    for j in 1:dims[2]
        for i in 1:dims[1]
            pint = potential_interactions[i,j,:]
            if sum(pint) > 0
                cell_poly = LG.Polygon(cell_coords(writer.xg[j], writer.xg[j+1], writer.yg[i], writer.yg[i+1]))
                cell_poly = if length(topography) > 0
                    topography_poly = LG.MultiPolygon(topography.coords)
                    LG.difference(cell_poly, topography_poly)
                end
                floeidx = collect(1:length(floes))[pint .== 1]
                # fic -> floes in cell - entire floe that is partially within grid cell
                # pic -> partially in cell - only includes pieces of floes that are within grid bounds
                pic_polys = [LG.intersection(cell_poly, LG.Polygon(floes[idx].coords)) for idx in floeidx]
                pic_area = [LG.area(poly) for poly in pic_polys]
                floeidx = floeidx[pic_area .> 0]
                pic_area = pic_area[pic_area .> 0]
                fic = floes[floeidx]

                floe_areas = fic.area
                area_ratios = pic_area ./ floe_areas
                area_tot = sum(pic_area)
                mass_tot = sum(fic.mass .* area_ratios)

                ds["mass"][istep, i, j] = mass_tot
                ds["area"][istep, i, j] = area_tot
                ds["si_frac"][istep, i, j] = area_tot/LG.area(cell_poly)
                # mass and area ratio
                ma_ratios = mass_tot > 0 ? area_ratios .* (fic.mass ./ mass_tot) : 0

                ds["u"][istep, i, j] = sum(fic.u .* ma_ratios)
                ds["v"][istep, i, j] = sum(fic.v .* ma_ratios)
                ds["du"][istep, i, j] = sum(fic.p_dudt .* ma_ratios)
                ds["dv"][istep, i, j] = sum(fic.p_dvdt .* ma_ratios)
                # need to add stress and strain!
                ds["overlap"][istep, i, j] = sum(fic.overarea)/length(fic)
                ds["height"][istep, i, j] = sum(fic.height .* ma_ratios)
            else
                ds["mass"][istep, i, j] = 0.0
                ds["area"][istep, i, j] = 0.0
                ds["si_frac"][istep, i, j] = 0.0
                ds["u"][istep, i, j] = 0.0
                ds["v"][istep, i, j] = 0.0
                ds["du"][istep, i, j] = 0.0
                ds["dv"][istep, i, j] = 0.0
                # need to add stress and strain!
                ds["overlap"][istep, i, j] = 0.0
                ds["height"][istep, i, j] = 0.0
            end
            
        end
    end
    close(ds)
end

"""
    setup_output_file!(writer::GridOutputWriter, nΔt, t::Type{T} = Float64)

Create output NetCDF file for given grid output writer. The file, which will have name defined in writer.fn will be saved in folder output/grid and will contain data for all fields defined by the writer. It will save data for nΔt/Δnout timesteps and the saved data will be of type T. 
Inputs:
        writer  <GridOutputWriter>
        nΔt     <Int> total number of timesteps in the simulation
        T       <Type> datatype to convert grid fields - must be a Float!
Outputs:
        Saved NetCDF file with needed fields and number of timesteps to be written to throughout the rest of the model. 
 """
function setup_output_file!(writer::GridOutputWriter, nΔt::Int, t::Type{T} = Float64) where T
    # Create file
    file_path = joinpath(pwd(), "output", "grid")
    !isdir(file_path) && mkdir(file_path)
    outfn = joinpath(file_path, writer.fn)
    isfile(outfn) && rm(outfn)
    ds = NCDataset(outfn, "c")
    ds.attrib["type"] = "Averaged grid data"

    # Define dimensions
    defDim(ds, "time", Inf)
    t = defVar(ds, "time", T, ("time",))
    t[:] = 0:writer.Δtout:nΔt
    t.attrib["units"] = "10 seconds"
    
    defDim(ds, "x", length(writer.xc))
    x = defVar(ds, "x", T, ("x",))
    x[:] = writer.xc
    defDim(ds, "y", length(writer.yc))
    y = defVar(ds, "y", T, ("y",))
    y[:] = writer.yc

    # Define variables
    for i in eachindex(writer.names)
        var = defVar(ds, writer.names[i], T, ("time", "y", "x"))
        var.attrib["units"] = writer.units[i]
        var.attrib["comments"] = writer.comments[i]
    end
    # Create file
    close(ds)
    return
end

"""
    setup_output_file!(writer::FloeOutputWriter, nΔt, t::Type{T} = Float64)

Create output NetCDF file for given floe output writer. The file, which will have name defined in writer.fn will be saved in folder output/floe and will contain data for all fields defined by the writer. It will save data for nΔt/Δnout timesteps and the saved data will be of type T. 
Inputs:
        writer  <FloeOutputWriter>
        nΔt     <Int> total number of timesteps in the simulation
        T       <Type> datatype to convert grid fields - must be a Float!
Outputs:
        Saved NetCDF file with needed fields and number of timesteps to be written to throughout the rest of the model. 
 """
function setup_output_file!(writer::FloeOutputWriter, nΔt, t::Type{T} = Float64) where T
    # Create file
    file_path = joinpath(pwd(), "output", "floe")
    !isdir(file_path) && mkdir(file_path)
    outfn = joinpath(file_path, writer.fn)
    isfile(outfn) && rm(outfn)
    ds = NCDataset(outfn, "c")
    ds.attrib["type"] = "Floe data"

    # Define dimensions
    defDim(ds, "time", Inf)
    t = defVar(ds, "time", T, ("time",))
    t.attrib["units"] = "10 seconds"
    t[:] = 0:writer.Δtout:nΔt

    defDim(ds, "floes", Inf)
    f = defVar(ds, "floes", T, ("floes",))
    f.attrib["units"] = "floeID"

    defDim(ds, "points", Inf)
    p = defVar(ds, "points", T, ("points",))
    p.attrib["units"] = "coordinate points"

    # Define variables
    for i in eachindex(writer.names)
        var = defVar(ds, writer.names[i], T, ("time", "floes"))
        var.attrib["units"] = writer.units[i]
        var.attrib["comments"] = writer.comments[i]
    end
    # Coords are seperate due to extra dimension
    xcoords = defVar(ds, "xcoords", T, ("time", "floes", "points"))
    xcoords.attrib["units"] = "location"
    xcoords.attrib["comments"] = "Floe x-coordinates"
    ycoords = defVar(ds, "ycoords", T, ("time", "floes", "points"))
    ycoords.attrib["units"] = "location"
    ycoords.attrib["comments"] = "Floe y-coordinates"

    # Close file
    close(ds)
    return
end

function write_data!(writer::GridOutputWriter, tstep, model)
    calc_eulerian_data!(model.floes, model.topos, writer, div(tstep, writer.Δtout) + 1)
end


function write_data!(writer::FloeOutputWriter, tstep, model)
    
    istep = div(tstep, writer.Δtout) + 1
    # I should make a Macro for this
    ds = NCDataset(joinpath(pwd(), "output", "floe", writer.fn), "a")

    nfloes = length(model.floes)
    Δfloes = nfloes - ds.dim["floes"]
    if Δfloes > 0
        ds["floes"][:] = 1:nfloes
        # might need to fill in previous values with missing or NaN??
    end

    ds["height"][istep, :] = model.floes.height
    ds["area"][istep, :] = model.floes.area
    ds["mass"][istep, :] = model.floes.mass
    ds["moment"][istep, :] = model.floes.moment
    ds["rmax"][istep, :] = model.floes.rmax
    ds["α"][istep, :] = model.floes.α
    ds["u"][istep, :] = model.floes.u
    ds["v"][istep, :] = model.floes.v
    ds["ξ"][istep, :] = model.floes.ξ
    ds["fxOA"][istep, :] = model.floes.fxOA
    ds["fyOA"][istep, :] = model.floes.fyOA
    ds["torqueOA"][istep, :] = model.floes.torqueOA
    ds["p_dxdt"][istep, :] = model.floes.p_dxdt
    ds["p_dydt"][istep, :] = model.floes.p_dydt
    ds["p_dudt"][istep, :] = model.floes.p_dudt
    ds["p_dvdt"][istep, :] = model.floes.p_dvdt
    ds["p_dξdt"][istep, :] = model.floes.p_dξdt
    ds["p_dαdt"][istep, :] = model.floes.p_dαdt
    ds["overarea"][istep, :] = model.floes.overarea
    ds["alive"][istep, :] = model.floes.alive
    xcentroid, ycentroid = seperate_xy([model.floes.centroid])
    ds["xcentroid"][istep, :] = xcentroid
    ds["ycentroid"][istep, :] = ycentroid
    for i in eachindex(model.floes)
        xcoords, ycoords = seperate_xy(model.floes[i].coords)
        ds["xcoords"][istep, i, 1:length(xcoords)] = xcoords
        ds["ycoords"][istep, i, 1:length(xcoords)] = ycoords
    end
    close(ds)
end

