"""
Structs and functions to calculate and write output from the simulation
"""

"""
    GridOutput

Options for averaged grid output. A list of these GridOutput objects is provided to the
GridOutputWriter constructor to determine which are saved during the model run.
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

#----------------------- Types of Output Writers -----------------------#
"""
    AbstractOutputWriter

An abstract type for output writers that provide data from simulation runs.
"""
abstract type AbstractOutputWriter end


"""
    CheckpointOutputWriter(Δtout, fn){ST<:AbstractString}<:AbstractOutputWriter

Checkpoint subtype of AbstractOutputWriter that holds information for outputting checkpoint information used for
restarting the simulation from a point where the writer saved data. Checkpoint data is saved every Δtout timesteps
to filepath. If the given file doesn't end in ".jld2" the extension will be appended. If overwrite_existing is true
then if there is already a file of the given name, it will be overwriten. Else it will thrown an error. 
"""
struct CheckpointOutputWriter<:AbstractOutputWriter
    Δtout::Int                  # Number of timesteps between checkpoint outputs
    filepath::String            # Filename for output file
    overwrite_existing::Bool    # Remove existing files if their filenames conflict.
    function CheckpointOutputWriter(Δtout, filepath, overwrite_existing)
        if isempty(filepath)
            throw(ArgumentError("A file path must be provided for each checkpoint writer."))
        end
        filepath = auto_extension(filepath, ".jld2")
        overwrite_existing && isfile(filepath) && rm(filepath, force=true)
        new(Δtout, filepath, overwrite_existing)
    end 
end

"""
    CheckpointOutputWriter(Δtout)

CheckpointOutputWriter writer that outputs need data to restart simulation at given timesteps Δtout
and saves the file in a file called "checkpoint.jld2".
Inputs:
        Δtout       <Int> number of timesteps between output
        dir         <String> Directory to save output to. Default: "." (current working directory).
        filename    <String> Descriptive filename.
        overwrtie
Outputs: FloeOutputWriter that outputs all Floe fields every Δtout timesteps to file.
"""
CheckpointOutputWriter(Δtout; dir = ".", filename = "checkpoint.jld2", overwrite = false) =
    CheckpointOutputWriter(Δtout, joinpath(dir, filename), overwrite)


"""
    FloeOutputWriter{ST<:AbstractString}<:AbstractOutputWriter

Floe subtype of AbstractOutputWriter that holds information for outputting floe information from model throughout simulation.
Output will be saved to the file defined by fn every Δtout timesteps.
Only outputs within the outputs list will be saved.
File will be saved as a JLD2 file to filepath. If the given file doesn't end in ".jld2" the extension will be appended.
If overwrite_existing is true then if there is already a file of the given name, it will be overwriten. Else it will thrown an error. 
"""
struct FloeOutputWriter<:AbstractOutputWriter
    outputs::Vector{Symbol}     # Floe fields to output
    Δtout::Int                  # Number of timesteps between floe outputs
    filepath::String            # Filename for output file
    overwrite_existing::Bool    # Remove existing files if their filenames conflict.

    function FloeOutputWriter(outputs, Δtout, filepath, overwrite_existing)
        if isempty(filepath)
            throw(ArgumentError("A filename must be provided for each floe writer."))
        end
        filepath = auto_extension(filepath, ".jld2")
        new(outputs, Δtout, filepath, overwrite_existing)
    end
end

"""
    FloeOutputWriter(Δtout, fn)

FloeOutput writer that outputs all Floe fields at given timesteps Δtout and saves the
information in a file of the provided name.
Inputs:
        Δtout   <Int> number of timesteps between output
        fn      <String> name of file to save grid data to
Outputs: FloeOutputWriter that outputs all Floe fields every Δtout timesteps to file fn
"""
FloeOutputWriter(Δtout; dir = ".", filename = "checkpoint.jld2", overwrite = false) =
    FloeOutputWriter(collect(fieldnames(Floe)), Δtout, joinpath(dir, filename), overwrite)


"""
    GridOutputWriter{FT<:AbstractFloat, ST<:AbstractString}<:AbstractOutputWriter

Grid subtype of AbstractOutputWriter that holds information for outputting Eularian data on the grid.
This output does not need to be the grid defined for the model. This grid can be coarser, or more fine, as defined by the
xg, yg, xc, and yc fields. Output on this scale will be saved to the file defined by fn every Δtout timesteps.
Data will be collected in the data field during calculation for easier writing to the NetCDF file.
Only outputs within the outputs list will be saved.
"""
struct GridOutputWriter{FT<:AbstractFloat, ST<:AbstractString}<:AbstractOutputWriter
    outputs::Vector{GridOutput}
    Δtout::Int              # Number of timesteps between grid outputs
    fn::ST                  # Filename for output file
    xg::Vector{FT}          # Grid lines in x-direction for ouput calculations
    yg::Vector{FT}          # Grid lines in y-direction for ouput calculations
    xc::Vector{FT}          # Center-lines on grid cells in the x-direction
    yc::Vector{FT}          # Center-lines on grid cells in the y-direction
    data::Array{FT, 3}      # Data to write to file stored here

    GridOutputWriter(outputs, Δtout, fn, xg, yg, xc, yc, data) = 
        length(xg) > 1 && length(yg) > 1 &&
        length(xg) - 1 == length(xc) == size(data)[2] &&
        length(yg) - 1 == length(yc) == size(data)[1] &&
        (length(outputs) == size(data)[3]) && length(fn) > 0 ?
        new{eltype(xg), typeof(fn)}(outputs, Δtout, fn, xg, yg, xc, yc, data) :
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
        t        <Type> datatype to convert saved data - must be a Float!
Output:
        GridOutputWriter that re-grids given grid to given dimensions, and number of timesteps between saving averaged output on this re-gridded grid to filename.
"""
function GridOutputWriter(outputs, Δtout, fn, grid::AbstractGrid, dims, ::Type{T} = Float64) where T
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
    data = zeros(T, length(yc), length(xc), length(outputs))
    return GridOutputWriter(outputs, Δtout, fn, xg, yg, xc, yc, data)
end

#----------------------- Initial Setup -----------------------#

""""
    output_initial_state!(sim::Simulation)

Write and save file that holds the initial state of the simulation.
Inputs:
        sim      <Simulation> simulation to run
        fn       <String> filename to save initial state file
Ouputs: File with initial simulation state will be saved in output/sim_name/initial_state.jld2.
"""
function output_initial_state!(sim::Simulation, fn = "initial_state.jld2")
    outfn = setup_directory!(fn, sim.name)
    jldsave(outfn; sim)
    return
end

"""
    setup_directory!(fn, sim_name)

Checks for a directory and creates it if not found.
Checks for file within director and deletes it if found. 
Inputs:
        fn          <String> name of file to check for
        sim_name    <String> name of simulation to be run
Outputs: Create directory output/sim_name if it doesn't exist and deletes the file output/sim_name/fn
if it exists. 
"""
function setup_directory!(fn, sim_name)
    # Create file and folder if needed
    file_path = joinpath("output", sim_name)
    !isdir(file_path) && mkdir(file_path)
    outfn = joinpath(file_path, fn)
    isfile(outfn) && rm(outfn)
    return outfn
end

"""
    setup_output_file!(writer::CheckpointOutputWriter, nΔt::Int, sim_name::String, ::Type{T} = Float64)

Sets up directory and file for given checkpoint writer to write to.
Inputs:
        writer      <CheckpointOutputWriter> writer to set up file for
        nΔt         <Int> number of timesteps to run - not used in this method 
        sim_name    <String> name of simulation
                    <Float> Float type simulation is run in - not used in this method
Outputs: Create directory output/sim_name if it doesn't exist and deletes the file output/sim_name/fn
if it exists. 
"""
function setup_output_file!(writer::CheckpointOutputWriter, nΔt::Int, sim_name::String, ::Type{T} = Float64) where T
    setup_directory!(writer.fn, sim_name)
    return
end

"""
    setup_output_file!(writer::FloeOutputWriter, nΔt::Int, sim_name::String, ::Type{T} = Float64)

Sets up directory and file for given floe writer to write to.
Inputs:
        writer      <FloeOutputWriter> writer to set up file for
        nΔt         <Int> number of timesteps to run - not used in this method 
        sim_name    <String> name of simulation
                    <Float> Float type simulation is run in - not used in this method
Outputs: Create directory output/sim_name if it doesn't exist and deletes the file output/sim_name/fn
if it exists. 
"""
function setup_output_file!(writer::FloeOutputWriter, nΔt::Int, sim_name::String, ::Type{T} = Float64) where T
    outfn = setup_directory!(writer.fn, sim_name)
    jldopen(outfn, "a+") do file
        for field in writer.outputs
            JLD2.Group(file, string(field))
        end
    end
end

"""
    setup_output_file!(writer::GridOutputWriter, nΔt::Int, sim_name::String, ::Type{T} = Float64) where T

Create output NetCDF file for given GridOutputWriter.
The file, named writer.fn, will be saved in folder output/sim_name and will contain data for fields defined by the writer.
It will save data for nΔt/Δnout timesteps and the saved data will be of type T. 
Inputs:
        writer      <GridOutputWriter>
        nΔt         <Int> total number of timesteps in the simulation
        sim_name    <String> simulation name to use for as folder name
        T           <Type> datatype to convert grid fields - must be a Float!
Outputs:
        Setup and saved NetCDF file with desired fields and number of timesteps to be written to throughout the rest of the simulation. 
 """
function setup_output_file!(writer::GridOutputWriter, nΔt::Int, sim_name::String, ::Type{T} = Float64) where T
    # Create file and folder if needed
    outfn = setup_directory!(writer.fn, sim_name)
    ds = NCDataset(outfn, "c")
    ds.attrib["type"] = "Gridded data"

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
        name = getname(output)
        unit, comment = getattrs(output)
        var = defVar(ds, name, T, ("time", "y", "x"))
        var.attrib["units"] = unit
        var.attrib["comments"] = comment
    end
    # Write to file and close
    close(ds)
    return
end

#----------------------- Write Data -----------------------#
"""
    write_data!(writer::CheckpointOutputWriter, tstep, model, sim_name)

Writes model's floe, ocean, and atmosphere data to JLD2 file. Data can be used to restart simulation run.

Inputs:
        writer      <FloeOutputWriter>
        tstep       <Int> simulation timestep
        model       <StructArray{Model}> model being simulated
        sim_name    <String> simulation name to use for as folder name
Output:
        Writes floes, ocean, and atmosphere to JLD2 file with name writer.fn for current timestep, which will be the group in the JLD2 file. 
"""
function write_data!(writer::CheckpointOutputWriter, tstep, model, sim_name)
    jldopen(joinpath("output", sim_name, writer.fn), "a+") do file
        timestep = JLD2.Group(file, string(tstep))
        timestep["floes"] = model.floes
        timestep["ocean"] = model.ocean
        timestep["atmosphere"] = model.atmos
    end
end

"""
    write_data!(writer::FloeOutputWriter, tstep, floes, sim_name)

Writes desired FloeOutputWriter data to JLD2 file.

Inputs:
        writer      <FloeOutputWriter>
        tstep       <Int> simulation timestep
        floes       <StructArray{Floe}> floes being simulated
        sim_name    <String> simulation name to use for as folder name
Output:
        Writes desired fields writer.outputs to JLD2 file with name writer.fn for current timestep, which will be the group in the JLD2 file. 
"""
function write_data!(writer::FloeOutputWriter, tstep, model, sim_name)
    jldopen(joinpath("output", sim_name, writer.fn), "a+") do file
        for output in writer.outputs
            file[string(output, "/", tstep)] = StructArrays.component(model.floes, output)
        end
    end
end

"""
    write_data!(writer::GridOutputWriter, tstep, model)

Writes GridOutputWriter data to NetCDF file created with setup_output_file! function.
Inputs:
        writer      <GridOutputWriter>
        tstep       <Int> simulation timestep
        model       <Model> model being simulated
        sim_name    <String> simulation name to use for as folder name
Output:
        Writes desired fields writer.outputs to file with name writer.fn for current timestep
"""
function write_data!(writer::GridOutputWriter, tstep, model, sim_name)
    live_floes = filter(f -> f.alive == 1, model.floes)
    if length(live_floes) > 0
        calc_eulerian_data!(live_floes, model.topos, writer)
        istep = div(tstep, writer.Δtout) + 1  # Julia indicies start at 1
        
        # Open file and write data from grid writer
        ds = NCDataset(joinpath(pwd(), "output", sim_name, writer.fn), "a+")
        for i in eachindex(writer.outputs)
            name = getname(writer.outputs[i])
            ds[name][istep, :, :] = writer.data[:, :, i]
        end
        close(ds)
    end
    return
end

#----------------------- Other Stuff -----------------------#
"""
    auto_extension(filename, ext) 

If `filename` ends in `ext`, return `filename`. Otherwise return `filename * ext`.
"""
function auto_extension(filename, ext) 
    Next = length(ext)
    filename[end-Next+1:end] == ext || (filename *= ext)
    return filename
end

"""
    rect_coords(xmin, xmax, ymin, ymax)

PolyVec coordinates of a rectangle given minimum and maximum x and y coordinates
Inputs:
    xmin    <Float> minimum x coordinate of rectangle
    xmax    <Float> maximum x coordinate of rectangle
    ymin    <Float> minimum y coordiante of rectangle
    ymax    <Float> maximum y coordiante of rectangle
Output:
    PolyVect coordinates for edges of rectangle with given minimums and maximums
"""
function rect_coords(xmin, xmax, ymin, ymax)
return [[[xmin, ymax], [xmin, ymin],
         [xmax, ymin], [xmax, ymax],
         [xmin, ymax]]]
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
        Floe data averaged on eularian grid provided and saved in writer.data field 
"""
function calc_eulerian_data!(live_floes, topography, writer)
    # Calculate/collect needed values
    Δx = writer.xg[2] - writer.xg[1]
    Δy = writer.yg[2] - writer.yg[1]
    cell_rmax = sqrt(Δx^2 + Δy^2)
    xgrid, ygrid = grids_from_lines(writer.xc, writer.yc)
    dims = size(xgrid)
    floe_centroids = live_floes.centroid
    floe_rmax = live_floes.rmax

    # Identify floes that potentially overlap each grid square and create mask
    potential_interactions = zeros(dims[1], dims[2], length(live_floes))
    for i in eachindex(live_floes)
        pint = sqrt.((xgrid .- floe_centroids[i][1]).^2 .+ (ygrid .- floe_centroids[i][2]).^2) .- (floe_rmax[i] + cell_rmax)
        pint[pint .> 0] .= 0
        pint[pint .< 0] .= 1
        potential_interactions[:,:,i] = pint
    end
    
    # Loop over each grid square
    for j in 1:dims[2]
        for i in 1:dims[1]
            pint = potential_interactions[i,j,:]
            # If there are any potential interactions
            if sum(pint) > 0
                cell_poly = LG.Polygon(rect_coords(writer.xg[j], writer.xg[j+1], writer.yg[i], writer.yg[i+1]))
                if length(topography) > 0
                    topography_poly = LG.MultiPolygon(topography.coords)
                    cell_poly = LG.difference(cell_poly, topography_poly)
                end

                floeidx = collect(1:length(live_floes))[pint .== 1]
                # fic -> floes in cell - entirety of floe that is partially within grid cell
                # pic -> partially in cell - only includes pieces of floes that are within grid bounds
                pic_polys = [LG.intersection(cell_poly, LG.Polygon(live_floes[idx].coords)) for idx in floeidx]
                pic_area = [LG.area(poly) for poly in pic_polys]
                floeidx = floeidx[pic_area .> 0]
                pic_area = pic_area[pic_area .> 0]
                fic = live_floes[floeidx]
                fic_area = fic.area

                area_ratios = pic_area ./ fic_area
                area_tot = sum(pic_area)
                mass_tot = sum(fic.mass .* area_ratios)

                if mass_tot > 0
                    # mass and area ratios
                    ma_ratios = area_ratios .* (fic.mass ./ mass_tot)
                    outputs = writer.outputs
                    for k in eachindex(outputs)
                        data = if outputs[k] == u_grid
                            sum(fic.u .* ma_ratios)
                        elseif outputs[k] == v_grid
                            sum(fic.v .* ma_ratios)
                        elseif outputs[k] == du_grid
                            sum(fic.p_dudt .* ma_ratios)
                        elseif outputs[k] == dv_grid
                            sum(fic.p_dvdt .* ma_ratios)
                        elseif outputs[k] == si_frac_grid
                            area_tot/LG.area(cell_poly)
                        elseif outputs[k] == overarea_grid
                            sum(fic.overarea)/length(fic)
                        elseif outputs[k] == mass_grid
                            mass_tot
                        elseif outputs[k] == area_grid
                            area_tot
                        elseif outputs[k] == height_grid
                            sum(fic.height .* ma_ratios)
                        end
                        # need to add stress and strain!
                        writer.data[i, j, k] = data
                    end
                else
                    writer.data[i, j, :] .= 0.0
                end
            else
                writer.data[i, j, :] .= 0.0
            end
        end
    end
    return
end

#----------------------- Metadata -----------------------#

"""
    getname(output::GridOutput)

Returns string name for given GridOutput object for writing to NetCDF
Input:
        output <GridOutput>
Output:
        <String> name for output object that is will be written as within output file
"""
function getname(output::GridOutput)
    n = Int(output)
    name = 
        n == 1 ? "u" :
        n == 2 ? "v" :
        n == 3 ? "du" :
        n == 4 ? "dv" :
        n == 5 ? "si_frac" :
        n == 6 ? "overlap" :
        n == 7 ? "mass" :
        n == 8 ? "area" :
        n == 9 ? "height" :
        throw(ArgumentError("Grid output provided is not known."))
    return name
end

"""
    getattrs(output::GridOutput)

Returns unit and comment attributes for each output type to be saved within output NetCDF file
Input:
        output<GridOutput>
Output:
        <Tuple(String, String)> tuple of string units and comments to be saved to output NetCDF file
"""
function getattrs(output::GridOutput)
    n = Int(output)
    unit, comment = 
        n == 1 ? ("m/s", "Average x-velocity of floes within grid cell") :
        n == 2 ? ("m/s", "Average y-velocity of floes within grid cell") :
        n == 3 ? ("m/s^2", "Average x-acceleration of floes within grid cell") :
        n == 4 ? ("m/s^2", "Average y-acceleration of floes within grid cell") :
        n == 5 ? ("unitless", "Fraction of grid cell covered by ice floes") :
        n == 6 ? ("m^2", "Average overlap area of floes per floe within grid cell") :
        n == 7 ? ("kg", "Total mass of floes within grid cell") :
        n == 8 ? ("m^2", "Total area of floes within grid cell") :
        n == 9 ? ("m", "Average height of floes within grid cell") :
        throw(ArgumentError("Grid output provided is not known."))
    return unit, comment
end

"""
    getattrs(output::FloeOutput)

Returns unit and comment attributes for each output type to be saved within output NetCDF file
Input:
        output<FloeOutput>
Output:
        <Tuple(String, String)> tuple of string units and comments to be saved to output NetCDF file
"""
function getattrs(output::Symbol)
    unit, comment =
        output == :centroid ? ("location", "Coordinates of centorid as [x y]") :
        output == :coords ? ("location", "Floe coordinates - either x or y") :
        output == :height ? ("m", "Floe height (uniform over floe)") :
        output == :area ? ("m^2", "Floe area") :
        output == :mass ? ("kg", "Floe mass") :
        output == :rmax ? ("m", "Maximum radius within floe") :
        output == :moment ? ("kg m^2", "Floe mass moment of intertia") :
        output == :angles ? ("kg m^2", "Floe interior angles") :
        output == :mc_x ? ("location", "x-coordinates for floe monte carlo points") :
        output == :mc_y ? ("location", "y-coordinates for floe monte carlo points") :
        output == :α ? ("radians", "Floe rotation since starting position") :
        output == :u ? ("m/s", "Floe x-direction velocity") :
        output == :v ? ("m/s", "Floe y-direction velocity") :
        output == :ξ ? ("rad/s", "Floe angular velocity") :
        output == :alive ? ("unitless", "Flag if floe is still active in simulation") :
        output == :id ? ("unitless", "Floe ID -  ghosts have negative ID of parents") :
        output == :ghosts ? ("unitless", "Index of floe's ghost floes in list") :
        output == :fxOA ? ("N", "X-directional forces on floe from the ocean and atmosphere") :
        output == :fyOA ? ("N", "Y-directional forces on floe from the ocean and atmosphere") :
        output == :trqOA ? ("N m", "Torque on floes from ocean and atmosphere") :
        output == :hflx ? ("W/m^2", "Average heatflux acting on the floe") : 
        output == :overarea ? ("m^2", "Overlap area of floe with other floes") :
        output == :collision_force ? ("N", "Total collision force from floes/topography/border in form [xf yf]") :
        output == :collision_torque ? ("N", "Total collision torque from floes/topography/border") :
        output == :interactions ? ("[unitless, N, N, location, location, N m, m^2]",
                                    "Matrix of floe's interactions with following columns: ID of floe interacted with,
                                    collision x-force on floe, collision y-force on floe, collision x-point, collision y-point,
                                    collision torque on floe, overlap with other floe") :
        output == :p_dxdt ? ("m/s", "Floe x-velocity from previous time step") :
        output == :p_dydt ? ("m/s", "Floe y-velocity from previous time step") :
        output == :p_dudt ? ("m/s^2", "Floe x-acceleration from previous time step") :
        output == :p_dvdt ? ("m/s^2", "Floe y-acceleration from previous time step") :
        output == :p_dξdt ? ("rad/s^2", "Floe angular acceleration from the previous timestep") :
        output == :p_dαdt ? ("rad/s", "Floe angular velocity from previous timestep") :
        ("", "") # if symbol isn't found, return empty attributes
    return unit, comment
end

# -----------Tools-------------------------------------------------------------
"""
    COULD BE TRANSFORMED INTO A TOOL TO CONVERT JLD2 to NetCDF

    setup_output_file!(writer::FloeOutputWriter, nΔt, t::Type{T} = Float64)

Create output NetCDF file for given floe output writer. The file, which will have name defined in writer.fn will be saved in folder output/floe and will contain data for all fields defined by the writer. It will save data for nΔt/Δnout timesteps and the saved data will be of type T. 
Inputs:
        writer      <FloeOutputWriter>
        nΔt         <Int> total number of timesteps in the simulation
        sim_name    <String> simulation name to use for as folder name
        T           <Type> datatype to convert grid fields - must be a Float!
Outputs:
        Saved NetCDF file with needed fields and number of timesteps to be written to throughout the rest of the model. 
 """
# function floes_to_netcdf(fn, ::Type{T} = Float64) where T
#     path, extension = split(fn, ".")
#     @assert extension == "jld2" "Can't run floes_to_netcdf on a non-jld2 file."
#     outfn = string(path, ".nc")
#     isfile(outfn) && rm(outfn)
#     # Load data and determine dimenstions
#     jld2_data = load(fn)
#     keynames = split.(keys(jld2_data), "/")
#     times = sort(parse.(Int, unique(first.(keynames))))
#     fields = unique(last.(keynames))

#     max_floes = 0
#     max_points = 0
#     for t in times
#         coords = jld2_data[string(t, "/coords")]
#         nfloes = length(coords)
#         if nfloes > max_floes
#             max_floes = nfloes
#         end
#         # we can parallelize this map with pmap?
#         npoints = maximum(map(c -> length(c[1]), coords))
#         if npoints > max_points
#             max_points = npoints
#         end
#     end

#     # Create NetCDF file
#     ds = NCDataset(outfn, "c")
#     ds.attrib["type"] = "Floe data"

#     # Define dimensions
#     defDim(ds, "time", Inf)
#     t = defVar(ds, "time", T, ("time",))
#     t.attrib["units"] = "10 seconds"
#     t[:] = times

#     defDim(ds, "floes", Inf)
#     f = defVar(ds, "floes", T, ("floes",))
#     t[:] = 1:max_floes

#     defDim(ds, "points", Inf)
#     p = defVar(ds, "points", T, ("points",))
#     p.attrib["units"] = "coordinate points"
#     t[:] = 1:npoints

#     # Define variables
#     for f in fields
#         unit, comment = getattrs(f)
#         if f == "coords"
#             var = defVar(ds, "xcoords", T, ("time", "floes", "points"))
#             var = defVar(ds, "ycoords", T, ("time", "floes", "points"))
#         else
#             var = defVar(ds, f, T, ("time", "floes"))
#         end
#         var.attrib["units"] = unit
#         var.attrib["comments"] = comment
#     end
# end

# function setup_output_file!(writer::FloeOutputWriter, nΔt, sim_name::String, ::Type{T} = Float64) where T
#     # Create file
#     file_path = joinpath(pwd(), "output", sim_name)
#     !isdir(file_path) && mkdir(file_path)
#     outfn = joinpath(file_path, writer.fn)
#     isfile(outfn) && rm(outfn)
#     ds = NCDataset(outfn, "c")
#     ds.attrib["type"] = "Floe data"

#     # Define dimensions
#     defDim(ds, "time", Inf)
#     t = defVar(ds, "time", T, ("time",))
#     t.attrib["units"] = "10 seconds"
#     t[:] = 0:writer.Δtout:nΔt

#     defDim(ds, "floes", Inf)
#     f = defVar(ds, "floes", T, ("floes",))
#     f.attrib["units"] = "floeID"

#     defDim(ds, "points", Inf)
#     p = defVar(ds, "points", T, ("points",))
#     p.attrib["units"] = "coordinate points"

#     # Define variables
#     for output in writer.outputs
#         name = getname(output)
#         unit, comment = getattrs(output)
#         if output != xcoords_floe && output != ycoords_floe
#             var = defVar(ds, name, T, ("time", "floes"))
#         else  # Coordinates have an extra dimension
#             var = defVar(ds, name, T, ("time", "floes", "points"))
#         end
#         var.attrib["units"] = unit
#         var.attrib["comments"] = comment
#     end
#     # Write to file and close
#     close(ds)
#     return
# end


"""
    write_data!(writer::FloeOutputWriter, tstep, model)

Writes FloeOutputWriter data to NetCDF file created with setup_output_file! function.
Inputs:
        writer  <FloeOutputWriter>
        tstep   <Int> simulation timestep
        model   <Model> model being simulated
        sim_name    <String> simulation name to use for as folder name
Output:
        Writes desired fields writer.outputs to file with name writer.fn for current timestep
"""
# function write_data!(writer::FloeOutputWriter, tstep, model, sim_name)
#     live_floes = filter(f -> f.alive == 1, model.floes)
#     nfloes = length(live_floes)
#     if nfloes > 0
#         istep = div(tstep, writer.Δtout) + 1  # Julia indicies start at 1
#         # Open file 
#         ds = NCDataset(joinpath(pwd(), "output", sim_name, writer.fn), "a+")

#         # Change in number of floes
#         Δfloes = nfloes - ds.dim["floes"]

#         # If less floes, add NaN values to make data square with previous data
#         nfloeNaN = Δfloes < 0 ? -Δfloes : 0

#         # If more floes, increase floes dimension
#         if Δfloes > 0
#             ds["floes"][:] = 1:nfloes
#         end

#         # Get data from all floes in Floe StructArray and transform if needed
#         for o in writer.outputs
#             data =
#                 if o == height_floe
#                     live_floes.height
#                 elseif o == area_floe
#                     live_floes.area
#                 elseif o == mass_floe
#                     live_floes.mass
#                 elseif o == moment_floe
#                     live_floes.moment
#                 elseif o == rmax_floe
#                     live_floes.rmax
#                 elseif o == α_floe
#                     live_floes.α
#                 elseif o == u_floe
#                     live_floes.u
#                 elseif o == v_floe
#                     live_floes.v
#                 elseif o == ξ_floe
#                     live_floes.ξ
#                 elseif o == fxOA_floe
#                     live_floes.fxOA
#                 elseif o == fyOA_floe
#                     live_floes.fyOA
#                 elseif o == trqOA_floe
#                     live_floes.trqOA
#                 elseif o == p_dxdt_floe
#                     live_floes.p_dxdt
#                 elseif o == p_dydt_floe
#                     live_floes.p_dydt
#                 elseif o == p_dudt_floe
#                     live_floes.p_dudt
#                 elseif o == p_dvdt_floe
#                     live_floes.p_dvdt
#                 elseif o == p_dξdt_floe
#                     live_floes.p_dξdt
#                 elseif o == p_dαdt_floe
#                     live_floes.p_dαdt
#                 elseif o == overarea_floe
#                     live_floes.overarea
#                 elseif o == alive_floe
#                     live_floes.alive
#                 elseif o == xcentroid_floe
#                     # Seperate x centroid data
#                     xcentroid, ycentroid = seperate_xy([live_floes.centroid])
#                     xcentroid
#                 elseif o == ycentroid_floe
#                     # Seperate y centroid data
#                     xcentroid, ycentroid = seperate_xy([live_floes.centroid])
#                     ycentroid
#                 elseif o == xcoords_floe
#                     # Seperate x-coordinate data and square coordinate length between floes by adding NaNs
#                     xcoords = [Subzero.seperate_xy(f.coords)[1] for f in live_floes]
#                     # If more points, increase point dimensions and square data
#                     npoints = [length(c) for c in xcoords]
#                     maxpoints = maximum(npoints)
#                     Δpoints = maxpoints - ds.dim["points"]
#                     if Δpoints > 0
#                         ds["points"][:] = 1:maxpoints
#                     end
#                     npointsNaN = [n < 0 ? -n : 0 for n in npoints .- ds.dim["points"]]
#                     square_xcoords = [vcat(xcoords[i], fill(NaN, npointsNaN[i])) for i in eachindex(xcoords)]
#                     # Rectangular matrix of floes by x-coordinate points
#                     hcat(square_xcoords...)'
#                 elseif o == ycoords_floe
#                     # Seperate y-coordinate data and square coordinate length between floes by adding NaNs
#                     ycoords = [Subzero.seperate_xy(f.coords)[2] for f in live_floes]
#                     # If more points, increase point dimensions and square data
#                     npoints = [length(c) for c in ycoords]
#                     maxpoints = maximum(npoints)
#                     Δpoints = maximum(npoints) - ds.dim["points"]
#                     if Δpoints > 0
#                         ds["points"][:] = 1:maxpoints
#                     end
#                     npointsNaN = [n < 0 ? -n : 0 for n in npoints .- ds.dim["points"]]
#                     square_ycoords = [vcat(ycoords[i], fill(NaN, npointsNaN[i])) for i in eachindex(ycoords)]
#                     # Rectangular matrix of floes by y-coordinate points
#                     hcat(square_ycoords...)'
#                 elseif o == xcollision_force
#                     first.(live_floes.collision_force)
#                 elseif o == ycollision_force
#                     last.(live_floes.collision_force)
#                 elseif o == collision_trq
#                     live_floes.collision_trq
#                 end
            
#             name = getname(o)
#             if length(size(data)) == 1  # if 1D (not coords)
#                 if Δfloes > 0
#                     ds[name][1:istep-1, (nfloes-Δfloes+1):nfloes] = fill(NaN, istep-1, Δfloes)
#                 end
#                 ds[name][istep, :] = vcat(data, fill(NaN, nfloeNaN))
#             elseif length(size(data)) == 2  # if 2D (coords)
#                 if Δfloes > 0
#                     ds[name][1:istep-1, (nfloes-Δfloes+1):nfloes, 1:ds.dim["points"]] =
#                         fill(NaN, istep-1, Δfloes, ds.dim["points"])
#                 end
#                 ds[name][istep, :, :] = vcat(data, fill(NaN, nfloeNaN, ds.dim["points"]))
#             end
            
#         end
#         close(ds)
#     end
#     return
# end