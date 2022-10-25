"""
Structs and functions to calculate and write output from the simulation
"""

@enum GridOutput begin
    u_grid = 1
    v_grid = 2
    du_grid = 3
    dv_grid = 4
    si_frac_grid = 5
    overarea_grid = 6
    mass_grid = 7
    area_grid = 8
    height_grid = 9
end

@enum FloeOutput begin
    height_floe = 1
    area_floe = 2
    mass_floe = 3
    moment_floe = 4
    rmax_floe = 5
    α_floe = 6
    xcentroid_floe = 7
    ycentroid_floe = 8
    u_floe = 9
    v_floe = 10
    ξ_floe = 11
    fxOA_floe = 12
    fyOA_floe = 13
    torqueOA_floe = 14
    p_dxdt_floe = 15
    p_dydt_floe = 16
    p_dudt_floe = 17
    p_dvdt_floe = 18
    p_dξdt_floe = 19
    p_dαdt_floe = 20
    overarea_floe = 21
    alive = 22
    xcoords = 23
    ycoords = 24
end

function getattrs(output::GridOutput)
    n = Int(output)
    name, unit, comment = 
        n == 1 ? ("u", "m/s", "Average x-velocity of floes within grid cell") :
        n == 2 ? ("v", "m/s", "Average y-velocity of floes within grid cell") :
        n == 3 ? ("du", "m/s^2", "Average x-acceleration of floes within grid cell") :
        n == 4 ? ("dv", "m/s^2", "Average y-acceleration of floes within grid cell") :
        n == 5 ? ("si_frac", "unitless", "Fraction of grid cell covered by ice floes") :
        n == 6 ? ("overlap", "m^2", "Average overlap area of floes per floe within grid cell") :
        n == 7 ? ("mass", "kg", "Total mass of floes within grid cell") :
        n == 8 ? ("area", "m^2", "Total area of floes within grid cell") :
        n == 9 ? ("height", "m", "Average height of floes within grid cell") :
        throw(ArgumentError("Grid output provided is not known."))
    return name, unit, comment
end

function getattrs(output::FloeOutput)
    n = Int(output)
    name, unit, comment =
        n == 1 ? ("height", "m", "Floe height (uniform over floe)") :
        n == 2 ? ("area", "m^2", "Floe area") :
        n == 3 ? ("mass", "kg", "Floe mass") :
        n == 4 ? ("moment", "kg m^2", "Floe mass moment of intertia") :
        n == 5 ? ("rmax", "m", "Maximum radius within floe") :
        n == 6 ? ("α", "radians", "Floe rotation since starting position") :
        n == 7 ? ("xcentroid", "location", "X-point of centorid") :
        n == 8 ? ("ycentroid", "location", "Y-point of centorid") :
        n == 9 ? ("u", "m/s", "Floe x-direction velocity") :
        n == 10 ? ("v", "m/s", "Floe y-direction velocity") :
        n == 11 ? ("ξ", "rad/s", "Floe angular velocity") :
        n == 12 ? ("fxOA", "N", "X-directional forces on floe from the ocean and atmosphere") :
        n == 13 ? ("fyOA", "N", "Y-directional forces on floe from the ocean and atmosphere") :
        n == 14 ? ("torqueOA", "N m", "Torque on floes from ocean and atmosphere") :
        n == 15 ? ("p_dxdt", "m/s", "Floe x-velocity from previous time step") :
        n == 16 ? ("p_dydt", "m/s", "Floe y-velocity from previous time step") :
        n == 17 ? ("p_dudt", "m/s^2", "Floe x-acceleration from previous time step") :
        n == 18 ? ("p_dvdt", "m/s^2", "Floe y-acceleration from previous time step") :
        n == 19 ? ("p_dξdt", "rad/s^2", "Floe angular acceleration from the previous timestep") :
        n == 20 ? ("p_dαdt", "rad/s", "Floe angular velocity from previous timestep") :
        n == 21 ? ("overarea", "m^2", "Overlap area of floe with other floes") :
        n == 22 ? ("alive", "unitless", "Flag if floe is still active in simulation") :
        n == 23 ? ("xcoords", "location", "Floe x-coordinates") :
        n == 24 ? ("ycoords", "location", "Floe y-coordinates") :
        throw(ArgumentError("Floe output provided is not known."))
    return name, unit, comment
end

"""
    GridOutputData{FT<:AbstractFloat}

Data averaged over a coarse grid covering the domain. Each slice of the array (i, j, slice) is the size of the output grid,
with one value per output grid square. This data is used for plotting and is saved to be able to examine model progress and data.
If all data is desired on the same resolution as the model, the coarse grid should be equal to the model grid and therefore the size
of each of these fields will be the same as the model grid's dimensions.
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

Grid subtype of AbstractOutputWriter that holds information for outputting Eularian data on the grid.
This output does not need to be the grid defined for the model. This grid can be coarser, or more fine as defined by the
xg, yg, xc, and yc fields. Output on this scale will be saved to the file defined by fn every Δtout timesteps.
"""
struct GridOutputWriter{FT<:AbstractFloat, ST<:AbstractString}<:AbstractOutputWriter
    outputs::Vector{GridOutput}
    Δtout::Int              # Number of timesteps between grid outputs
    fn::ST                  # Filename for output file
    xg::Vector{FT}          # Grid lines in x-direction for ouput calculations
    yg::Vector{FT}          # Grid lines in y-direction for ouput calculations
    xc::Vector{FT}          # Center-lines on grid cells in the x-direction
    yc::Vector{FT}          # Center-lines on grid cells in the y-direction
    data_arr::GridOutputData{FT}

    GridOutputWriter(outputs, Δtout, fn, xg, yg, xc, yc, data_arr) = 
        length(xg) > 1 && length(yg) > 1 &&
        length(xg) - 1 == length(xc) == size(data_arr.data)[2] &&
        length(yg) - 1 == length(yc) == size(data_arr.data)[1] &&
        (length(outputs) == size(data_arr.data)[3]) && length(fn) > 0 ?
        new{eltype(xg), typeof(fn)}(outputs, Δtout, fn, xg, yg, xc, yc, data_arr) :
        throw(ArgumentError("Output grid lines must have at least one grid square and number of outputs must match data array dimensions"))
end


"""
    GridOutputWriter(outputs, Δtout, fn, grid::Grid, dims)

Create GridOutputWriter for given grid, desired output dimensions to re-grid given grid, and desired number of timesteps, saved to given filename.
Inputs:
        outputs  <GridOutputs[]> list of GridOutputs desired for output
        Δtout    <Int> number of timesteps between output
        fn       <String> name of file to save grid data to
        grid     <Grid> original grid, which we are re-gridding
        dims     <(Int, Int)> output grid dimensions - rows -> ny, cols -> nx
        T        <Type> datatype to convert saved data - must be a Float!
Output:
        GridOutputWriter that re-grids given grid to given dimensions, and number of timesteps between saving averaged output on this re-gridded grid to filename.
"""
function GridOutputWriter(outputs, Δtout, fn, grid::Grid, dims, t::Type{T} = Float64) where T
    # Define output grid
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
    # Data output container
    data = GridOutputData(zeros(T, length(yc), length(xc), length(outputs)))
    return GridOutputWriter(outputs, Δtout, fn, xg, yg, xc, yc, data)
end

"""
    FloeOutputWriter{FT<:AbstractFloat, ST<:AbstractString}<:AbstractOutputWriter

Floe subtype of AbstractOutputWriter that holds information for outputting floe information from model throughout simulation.
The edges of the grid are defined by xg and yg (the grid lines for the model) for plotting and reconstruction purposes. 
Output on this scale will be saved to the file defined by fn every Δtout timesteps.
"""
struct FloeOutputWriter{FT<:AbstractFloat, ST<:AbstractString}<:AbstractOutputWriter
    outputs::Vector{FloeOutput}
    Δtout::Int              # Number of timesteps between grid outputs
    fn::ST                  # Filename for output file
    xg::Vector{FT}          # Grid lines in x-direction for plotting
    yg::Vector{FT}          # Grid lines in y-direction for plotting

    FloeOutputWriter(outputs, Δtout, fn, xg, yg) = 
        length(xg) > 1 && length(yg) > 1 && length(fn) > 0  ?
        new{eltype(xg), typeof(fn)}(outputs, Δtout, fn, xg, yg) :
        throw(ArgumentError("Output grid lines must have at least one grid square and a filename must be provided."))
end

"""
    FloeOutputWriter(outputs, Δtout, fn, grid::Grid)

Create FloeOutputWriter for given grid and desired number of timesteps, saved to the given file name.
Inputs:
    outputs  <Vector{FloeOutputs}>
    Δtout    <Int> number of timesteps between output
    fn       <String> name of file to save grid data to
    grid     <Grid> original grid, which we are re-gridding
Output:
    GridOutputWriter that holds desired floe fields to save to given filename and the number of timesteps between saving this output.
    Grid extents is saved as metadata for plotting.
"""
FloeOutputWriter(outputs, Δtout, fn, grid::Grid) = FloeOutputWriter(outputs, Δtout, fn, grid.xg, grid.yg)

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
    for output in writer.outputs
        name, unit, comment = getattrs(output)
        var = defVar(ds, name, T, ("time", "y", "x"))
        var.attrib["units"] = unit
        var.attrib["comments"] = comment
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
    for output in writer.outputs
        name, unit, comment = getattrs(output)
        if output != xcoords && output != ycoords
            var = defVar(ds, name, T, ("time", "floes"))
        else
            var = defVar(ds, name, T, ("time", "floes", "points"))
        end
        var.attrib["units"] = unit
        var.attrib["comments"] = comment
    end
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

                    writer.data_arr.data[i, j, Int(u_grid)] = sum(fic.u .* ma_ratios)
                    writer.data_arr.data[i, j, Int(v_grid)] = sum(fic.v .* ma_ratios)
                    writer.data_arr.data[i, j, Int(du_grid)] = sum(fic.p_dudt .* ma_ratios)
                    writer.data_arr.data[i, j, Int(dv_grid)] = sum(fic.p_dvdt .* ma_ratios)
                    writer.data_arr.data[i, j, Int(si_frac_grid)] = area_tot/LG.area(cell_poly)
                    writer.data_arr.data[i, j, Int(overarea_grid)] = sum(fic.overarea)/length(fic)
                    writer.data_arr.data[i, j, Int(mass_grid)] = mass_tot
                    writer.data_arr.data[i, j, Int(area_grid)] = area_tot
                    writer.data_arr.data[i, j, Int(height_grid)] = sum(fic.height .* ma_ratios)
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
    ds["u"][istep, :, :] = writer.data_arr.data[:, :, Int(u_grid)]
    ds["v"][istep, :, :] = writer.data_arr.data[:, :, Int(v_grid)]
    ds["du"][istep, :, :] = writer.data_arr.data[:, :, Int(du_grid)]
    ds["dv"][istep, :, :] = writer.data_arr.data[:, :, Int(dv_grid)]
    ds["si_frac"][istep, :, :] = writer.data_arr.data[:, :, Int(si_frac_grid)]
    ds["overlap"][istep, :, :] = writer.data_arr.data[:, :, Int(overarea_grid)]
    ds["mass"][istep, :, :] = writer.data_arr.data[:, :, Int(mass_grid)]
    ds["area"][istep, :, :] = writer.data_arr.data[:, :, Int(area_grid)]
    ds["height"][istep, :, :] = writer.data_arr.data[:, :, Int(height_grid)]
    # need to add stress and strain!
    close(ds)
end


function write_data!(writer::FloeOutputWriter, tstep, model)
    istep = div(tstep, writer.Δtout) + 1  # Julia indicies start at 1
    # Open file and write data from floes
    ds = NCDataset(joinpath(pwd(), "output", "floe", writer.fn), "a")
    # Change in number of floes
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

