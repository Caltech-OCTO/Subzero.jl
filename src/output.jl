"""
Structs and functions to calculate and write output from the simulation
"""

"""
    GridOutputData{FT<:AbstractFloat}

Data averaged over a coarse grid covering the domain. Each slice of the array (i, j, slice) is the size of the output grid, with one value per output grid square. This data is used for plotting and is saved to be able to examine model progress and data. If all data is desired on the same resolution as the model, the coarse grid should be equal to the model grid and therefore the size of each of these fields will be the same as the model grid's dimensions.
"""
struct GridOutputData{FT<:AbstractFloat}
    data::Array{FT, 3}
    # Can add average flag and average array
end

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
    data_arr::GridOutputData{FT}

    GridOutputWriter(Δtout, fn, xg, yg, xc, yc, names, units, comments, data_arr) = 
        length(xg) > 1 && length(yg) > 1 &&
        length(xg) - 1 == length(xc) == size(data_arr.data)[2] &&
        length(yg) - 1 == length(yc) == size(data_arr.data)[1] &&
        (length(names) == length(units) == length(comments) == size(data_arr.data)[3]) && length(fn) > 0 ?
        new{eltype(xg), eltype(names)}(Δtout, fn, xg, yg, xc, yc, names, units, comments, data_arr) :
        throw(ArgumentError("Output grid lines must have at least one grid square or names, units, and comments must be the same size."))
end


"""
    GridOutputWriter(Δtout, fn, grid::Grid, dims)

Create GridOutputWriter for given grid, desired output dimensions to re-grid given grid, and desired number of timesteps, saved to given filename.
Inputs:
        Δtout    <Int> number of timesteps between output
        fn       <String> name of file to save grid data to
        grid     <Grid> original grid, which we are re-gridding
        dims     <(Int, Int)> output grid dimensions - rows -> ny, cols -> nx
        T        <Type> datatype to convert saved data - must be a Float!
Output:
        GridOutputWriter that re-grids given grid to given dimensions, and number of timesteps between saving averaged output on this re-gridded grid to filename.
"""
function GridOutputWriter(Δtout, fn, grid::Grid, dims, t::Type{T} = Float64) where T
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
    data = GridOutputData(zeros(T, length(yc), length(xc), length(grid_names)))
    return GridOutputWriter(Δtout, fn, xg, yg, xc, yc, grid_names, grid_units, grid_comments, data)
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
        length(xg) > 1 && length(yg) > 1 && (length(names) == length(units) == length(comments)) && length(fn) > 0  ?
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
function calc_eulerian_data!(floes, topography, writer)
    # Calculate needed values
    live_floes = filter(f -> f.alive == 1, floes)
    Δx = writer.xg[2] - writer.xg[1]
    Δy = writer.yg[2] - writer.yg[1]
    cell_rmax = sqrt(Δx^2 + Δy^2)
    xgrid, ygrid = grids_from_lines(writer.xc, writer.yc)
    dims = size(xgrid)
    floe_centroids = live_floes.centroid
    floe_rmax = live_floes.rmax

    # Identify floes that potentially overlap each grid square
    potential_interactions = zeros(dims[1], dims[2], length(floes))
    for i in eachindex(floes)
        pint = sqrt.((xgrid .- floe_centroids[i][1]).^2 .+ (ygrid .-
                      floe_centroids[i][2]).^2) .- (floe_rmax[i] + cell_rmax)
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

                if mass_tot > 0
                    # mass and area ratios
                    ma_ratios = area_ratios .* (fic.mass ./ mass_tot)

                    writer.data_arr.data[i, j, 1] = sum(fic.u .* ma_ratios)  # u
                    writer.data_arr.data[i, j, 2] = sum(fic.v .* ma_ratios)  # v
                    writer.data_arr.data[i, j, 3] = sum(fic.p_dudt .* ma_ratios)  # du
                    writer.data_arr.data[i, j, 4] = sum(fic.p_dvdt .* ma_ratios)  # dv
                    writer.data_arr.data[i, j, 5] = area_tot/LG.area(cell_poly)  # si_f
                    writer.data_arr.data[i, j, 6] = sum(fic.overarea)/length(fic)  #over
                    writer.data_arr.data[i, j, 7] = mass_tot  # mass
                    writer.data_arr.data[i, j, 8] = area_tot  # area
                    writer.data_arr.data[i, j, 9] = sum(fic.height .* ma_ratios) #height
                    # need to add stress and strain!
                else
                    fill!(writer.data_arr.data[i, j, :], 0.0)
                end
            else
                fill!(writer.data_arr.data[i, j, :], 0.0)
            end
        end
    end
end


function write_data!(writer::GridOutputWriter, tstep, model)
    calc_eulerian_data!(model.floes, model.topos, writer)

    istep = div(tstep, writer.Δtout) + 1  # Julia indicies start at 1
    # Open file and write data from grid writer
    ds = NCDataset(joinpath(pwd(), "output", "grid", writer.fn), "a")
    ds["u"][istep, :, :] = writer.data_arr.data[:, :, 1]
    ds["v"][istep, :, :] = writer.data_arr.data[:, :, 2]
    ds["du"][istep, :, :] = writer.data_arr.data[:, :, 3]
    ds["dv"][istep, :, :] = writer.data_arr.data[:, :, 4]
    ds["si_frac"][istep, :, :] = writer.data_arr.data[:, :, 5]
    ds["overlap"][istep, :, :] = writer.data_arr.data[:, :, 6]
    ds["mass"][istep, :, :] = writer.data_arr.data[:, :, 7]
    ds["area"][istep, :, :] = writer.data_arr.data[:, :, 8]
    ds["height"][istep, :, :] = writer.data_arr.data[:, :, 9]
    # need to add stress and strain!
    close(ds)
end


function write_data!(writer::FloeOutputWriter, tstep, model)
    
    istep = div(tstep, writer.Δtout) + 1  # Julia indicies start at 1
    # I should make a Macro for this
    ds = NCDataset(joinpath(pwd(), "output", "floe", writer.fn), "a")

    nfloes = length(model.floes)
    Δfloes = nfloes - ds.dim["floes"]
    if Δfloes > 0
        ds["floes"][:] = 1:nfloes
        # might need to fill in previous values with missing or NaN??
    end
    nNaN = Δfloes < 0 ? -Δfloes : 0

    ds["height"][istep, :] = vcat(model.floes.height, fill(NaN, nNaN))
    ds["area"][istep, :] = vcat(model.floes.area, fill(NaN, nNaN))
    ds["mass"][istep, :] = vcat(model.floes.mass, fill(NaN, nNaN))
    ds["moment"][istep, :] = vcat(model.floes.moment, fill(NaN, nNaN))
    ds["rmax"][istep, :] = vcat(model.floes.rmax, fill(NaN, nNaN))
    ds["α"][istep, :] = vcat(model.floes.α, fill(NaN, nNaN))
    ds["u"][istep, :] = vcat(model.floes.u, fill(NaN, nNaN))
    ds["v"][istep, :] = vcat(model.floes.v, fill(NaN, nNaN))
    ds["ξ"][istep, :] = vcat(model.floes.ξ, fill(NaN, nNaN))
    ds["fxOA"][istep, :] = vcat(model.floes.fxOA, fill(NaN, nNaN))
    ds["fyOA"][istep, :] = vcat(model.floes.fyOA, fill(NaN, nNaN))
    ds["torqueOA"][istep, :] = vcat(model.floes.torqueOA, fill(NaN, nNaN))
    ds["p_dxdt"][istep, :] = vcat(model.floes.p_dxdt, fill(NaN, nNaN))
    ds["p_dydt"][istep, :] = vcat(model.floes.p_dydt, fill(NaN, nNaN))
    ds["p_dudt"][istep, :] = vcat(model.floes.p_dudt, fill(NaN, nNaN))
    ds["p_dvdt"][istep, :] = vcat(model.floes.p_dvdt, fill(NaN, nNaN))
    ds["p_dξdt"][istep, :] = vcat(model.floes.p_dξdt, fill(NaN, nNaN))
    ds["p_dαdt"][istep, :] = vcat(model.floes.p_dαdt, fill(NaN, nNaN))
    ds["overarea"][istep, :] = vcat(model.floes.overarea, fill(NaN, nNaN))
    ds["alive"][istep, :] = vcat(model.floes.alive, fill(NaN, nNaN))
    xcentroid, ycentroid = seperate_xy([model.floes.centroid])
    ds["xcentroid"][istep, :] = vcat(xcentroid, fill(NaN, nNaN))
    ds["ycentroid"][istep, :] = vcat(ycentroid, fill(NaN, nNaN))

    for i in eachindex(model.floes)
        xcoords, ycoords = seperate_xy(model.floes[i].coords)
        npoints = length(xcoords)
        Δpoints = npoints - ds.dim["points"]
        if Δpoints > 0
            ds["points"][:] = 1:npoints
        end
        npNaN = Δpoints < 0 ? -Δpoints : 0
        ds["xcoords"][istep, i, :] = vcat(xcoords, fill(NaN, npNaN))
        ds["ycoords"][istep, i, :] = vcat(ycoords, fill(NaN, npNaN))
    end
    close(ds)
end

