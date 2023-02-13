"""
Structs and functions to calculate and write output from the simulation
"""

#----------------------- Types of Output Writers -----------------------#
"""
    AbstractOutputWriter

An abstract type for output writers that provide data from simulation runs.
"""
abstract type AbstractOutputWriter end

"""
    InitialStateOuputWriter<:AbstractOutputWriter
"""
struct InitialStateOutputWriter<:AbstractOutputWriter
    filepath::String            # Filename for output file
    overwrite::Bool             # Remove existing files if their filenames conflict.
end

function InitialStateOutputWriter(;dir = ".", filename = "initial_state.jld2", overwrite = false, jld2_kw = Dict{Symbol, Any}())
    filepath = initialize_jld2_file!(dir, filename, overwrite, [], jld2_kw)
    return InitialStateOutputWriter(filepath, overwrite)
end

"""
    CheckpointOutputWriter(Δtout, fn){ST<:AbstractString}<:AbstractOutputWriter

Checkpoint subtype of AbstractOutputWriter that holds information for outputting checkpoint information used for
restarting the simulation from a point where the writer saved data. Checkpoint data is saved every Δtout timesteps
to filepath. If the given file doesn't end in ".jld2" the extension will be appended. If overwrite is true
then if there is already a file of the given name, it will be overwriten. Else it will thrown an error. 
"""
struct CheckpointOutputWriter<:AbstractOutputWriter
    Δtout::Int                  # Number of timesteps between checkpoint outputs
    filepath::String            # Filename for output file
    overwrite::Bool             # Remove existing files if their filenames conflict.
end

"""
    CheckpointOutputWriter(Δtout)

CheckpointOutputWriter writer that outputs need data to restart simulation at given timesteps Δtout
and saves the file in a file called "checkpoint.jld2".
Inputs:
        Δtout       <Int> number of timesteps between output
        dir         <String> Directory to save output to. Default: "." (current working directory).
        filename    <String> Descriptive filename.
        overwrite
Outputs: FloeOutputWriter that outputs all Floe fields every Δtout timesteps to file.
"""
function CheckpointOutputWriter(Δtout; dir = ".", filename = "checkpoint.jld2", overwrite = false, jld2_kw = Dict{Symbol, Any}())
    filepath = initialize_jld2_file!(dir, filename, overwrite, [:ocean, :atmos, :floes], jld2_kw)
    return CheckpointOutputWriter(Δtout, filepath, overwrite)
end

"""
    FloeOutputWriter{ST<:AbstractString}<:AbstractOutputWriter

Floe subtype of AbstractOutputWriter that holds information for outputting floe information from model throughout simulation.
Output will be saved to the file defined by fn every Δtout timesteps.
Only outputs within the outputs list will be saved.
File will be saved as a JLD2 file to filepath. If the given file doesn't end in ".jld2" the extension will be appended.
If overwrite is true then if there is already a file of the given name, it will be overwriten. Else it will thrown an error. 
"""
struct FloeOutputWriter<:AbstractOutputWriter
    outputs::Vector{Symbol}     # Floe fields to output
    Δtout::Int                  # Number of timesteps between floe outputs
    filepath::String            # Filename for output file
    overwrite::Bool    # Remove existing files if their filenames conflict
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
function FloeOutputWriter(outputs, Δtout; dir = ".", filename = "floes.jld2", overwrite = false, jld2_kw = Dict{Symbol, Any}())
    filepath = initialize_jld2_file!(dir, filename, overwrite, outputs, jld2_kw)
    return FloeOutputWriter(outputs, Δtout, filepath, overwrite)
end

FloeOutputWriter(Δtout; dir = ".", filename = "floes.jld2", overwrite = false, jld2_kw = Dict{Symbol, Any}()) = 
    FloeOutputWriter(collect(fieldnames(Floe)), Δtout; dir = dir, filename = filename, overwrite = overwrite, jld2_kw = jld2_kw)


"""
    GridOutputWriter{FT<:AbstractFloat}<:AbstractOutputWriter

Grid subtype of AbstractOutputWriter that holds information for outputting floe data on the grid.
This output does not need to be the grid defined for the model. This grid can be coarser, or more fine, as defined by the
xg and yg fields. Output on this scale will be saved to the file defined by filepath every Δtout timesteps.
Data will be collected in the data field during calculation for easier writing to the NetCDF file.
Only outputs within the outputs list will be saved. There is a limited number of
floe outputs that can be calculated by the GridOutputWriter.
"""
struct GridOutputWriter{FT<:AbstractFloat}<:AbstractOutputWriter
    outputs::Vector{Symbol}
    Δtout::Int                  # Number of timesteps between grid outputs
    filepath::String                # Filename for output file
    overwrite::Bool    # Remove existing files if their filenames conflict
    xg::Vector{FT}              # Grid lines in x-direction for ouput calculations
    yg::Vector{FT}              # Grid lines in y-direction for ouput calculations
    data::Array{FT, 3}          # Data to write to file stored here
    average::Bool
end

function get_known_grid_outputs()
    return Set([:u_grid, :v_grid, :dudt_grid, :dvdt_grid, :overarea_grid, :mass_grid,
                :area_grid, :height_grid, :si_frac_grid, :stress_xx_grid, :stress_yx_grid, 
                :stress_xy_grid, :stress_yy_grid, :stress_eig_grid, :strain_ux_grid, :strain_vx_grid, :strain_uy_grid,
                :strain_vy_grid])
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
function GridOutputWriter(outputs::Vector{Symbol}, Δtout, grid::AbstractGrid, dims; dir = ".", filename = "gridded_data.nc", overwrite = false, average = false, t::Type{T} = Float64) where T
    # Check for known outputs - need routine to calculate in calc_eulerian_data
    known_grid_outputs = get_known_grid_outputs()
    remove_idx = []
    for i in eachindex(outputs)
        if !(outputs[i] in known_grid_outputs)
            @warn "$(outputs[i]) is not a known grid output. It will not be output with the GridOutputWriter."
            push!(remove_idx, i)
        end
    end
    deleteat!(outputs, remove_idx)

    # Define output grid
    lx = grid.xg[1]
    ux = grid.xg[end]
    ly = grid.yg[1]
    uy = grid.yg[end]
    xg = collect(lx:(ux-lx)/dims[2]:ux)
    yg = collect(ly:(uy-ly)/dims[1]:uy)

    # Create file path
    filepath = initialize_netcdf_file!(dir, filename, overwrite, outputs, xg, yg, T)

    # Data output container
    data = zeros(T, length(yg) - 1, length(xg) - 1, length(outputs))
    return GridOutputWriter(outputs, Δtout, filepath, overwrite, xg, yg, data, average)
end

GridOutputWriter(Δtout::Int, grid::AbstractGrid, dims; dir = ".", filename = "gridded_data.nc", overwrite = false, average = false, t::Type{T} = Float64) where T =
    GridOutputWriter(collect(get_known_grid_outputs()), Δtout, grid, dims, dir = dir,
                    filename = filename, overwrite = overwrite, average = average, t = T)

#----------------------- Write Data -----------------------#

function write_data!(writer::InitialStateOutputWriter, tstep, sim)
    jldopen(writer.filepath, "a") do file
        file["sim"] = sim
    end
end

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
function write_data!(writer::CheckpointOutputWriter, tstep, sim)
    jldopen(writer.filepath, "a+") do file
        file[string("floes/", tstep)] = sim.model.floes
        file[string("ocean/", tstep)] = sim.model.ocean
        file[string("atmos/", tstep)] = sim.model.atmos
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
function write_data!(writer::FloeOutputWriter, tstep, sim)
    jldopen(writer.filepath, "a+") do file
        for output in writer.outputs
            file[string(output, "/", tstep)] = StructArrays.component(sim.model.floes, output)
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
function write_data!(writer::GridOutputWriter, tstep, sim)
    live_floes = filter(f -> f.alive, sim.model.floes)
    if length(live_floes) > 0
        calc_eulerian_data!(live_floes, sim.model.domain.topography, writer)
        istep = div(tstep, writer.Δtout) + 1  # Julia indicies start at 1
        
        # Open file and write data from grid writer
        ds = NCDataset(writer.filepath, "a")
        for i in eachindex(writer.outputs)
            name = string(writer.outputs[i])
            ds[name][istep, :, :] = writer.data[:, :, i]
        end
        close(ds)
    end
    return
end

#----------------------- File Setup -----------------------#
"""Checks for a directory and creates it if not found.
    Checks for file within director and deletes it if found. 
    Inputs:
            fn          <String> name of file to check for
            sim_name    <String> name of simulation to be run
    Outputs: Create directory output/sim_name if it doesn't exist and deletes the file output/sim_name/fn
    if it exists. 
    """
function initialize_jld2_file!(dir, filename, overwrite, outputs, jld2_kw)
    # create path to file
    mkpath(dir)
    filename = auto_extension(filename, ".jld2")
    filepath = joinpath(dir, filename)
    # check if file already exists
    if !overwrite && isfile(filepath)
        throw(ErrorException("File $filepath already exists and overwrite is false."))
    end
    overwrite && isfile(filepath) && rm(filepath, force=true)
    try
        # Create the file and the needed groups
        jldopen(filepath, "a+"; jld2_kw...) do file
            # Add groups for each output and each output's metadata
            for o in outputs
                JLD2.Group(file, string(o))
                # Add metadata
                attrs = getattrs(o)
                file[string("metadata/", o)] = attrs
            end
        end
    catch err
        @warn """Initialization of $filepath failed because of $(sprint(showerror, err))"""
    end
    return filepath
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
function initialize_netcdf_file!(dir, filename, overwrite, outputs, xg, yg, T)
    mkpath(dir)
    filename = auto_extension(filename, ".nc")
    filepath = joinpath(dir, filename)
    if !overwrite && isfile(filepath)
        throw(ErrorException("File $filepath already exists and overwrite is false."))
    end
    overwrite && isfile(filepath) && rm(filepath, force=true)
    try
        # Create the file and the needed groups
        dataset = NCDataset(filepath, "c")
        dataset.attrib["type"] = "Floe data averaged on the grid. The grid is broken down into user provided dimensions."

        # Define dimensions
        defDim(dataset, "time", Inf)
        defVar(dataset, "time", T, ("time",), attrib = Dict("units" => "10 seconds"))

        defDim(dataset, "x", length(xg) - 1)
        x = defVar(dataset, "x", T, ("x",), attrib = Dict("units" => "meters"))
        x[:] = xg[1:end-1] .+ 0.5(xg[2] - xg[1])

        defDim(dataset, "y", length(yg) - 1)
        y = defVar(dataset, "y", T, ("y",), attrib = Dict("units" => "meters"))
        y[:] = yg[1:end-1] .+ 0.5(yg[2] - yg[1])

        # Define variables
        for o in outputs
            unit, comment = getattrs(o)
            defVar(dataset, string(o), T, ("time", "y", "x"), attrib = Dict("units" => unit, "comments" => comment))
        end
        # Write to file and close
        close(dataset)
    catch err
        @warn """Initialization of $filepath failed because of $(sprint(showerror, err))"""
    end
    return filepath
end

#----------------------- Utils -----------------------#
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
return [[[xmin, ymin], [xmin, ymax],
         [xmax, ymax], [xmax, ymin],
         [xmin, ymin]]]
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
function calc_eulerian_data!(floes, topography, writer)
    # Calculate/collect needed values
    Δx = writer.xg[2] - writer.xg[1]
    Δy = writer.yg[2] - writer.yg[1]
    cell_rmax = sqrt(Δx^2 + Δy^2)
    xgrid, ygrid = grids_from_lines(writer.xg[1:end-1] .+ 0.5Δx, writer.yg[1:end-1] .+ 0.5Δy)
    dims = size(xgrid)
    floe_centroids = floes.centroid
    floe_rmax = floes.rmax

    # Identify floes that potentially overlap each grid square and create mask
    potential_interactions = zeros(dims[1], dims[2], length(floes))
    for i in eachindex(floes)
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

                floeidx = collect(1:length(floes))[pint .== 1]
                # fic -> floes in cell - entirety of floe that is partially within grid cell
                # pic -> partially in cell - only includes pieces of floes that are within grid bounds
                pic_polys = [LG.intersection(cell_poly, LG.Polygon(floes[idx].coords)) for idx in floeidx]
                pic_area = [LG.area(poly) for poly in pic_polys]
                floeidx = floeidx[pic_area .> 0]
                pic_area = pic_area[pic_area .> 0]
                fic = floes[floeidx]
                fic_area = fic.area

                area_ratios = pic_area ./ fic_area
                area_tot = sum(pic_area)
                mass_tot = sum(fic.mass .* area_ratios)

                if mass_tot > 0
                    # mass and area ratios
                    ma_ratios = area_ratios .* (fic.mass ./ mass_tot)
                    outputs = writer.outputs
                    for k in eachindex(outputs)
                        data = if outputs[k] == :u_grid
                            sum(fic.u .* ma_ratios)
                        elseif outputs[k] == :v_grid
                            sum(fic.v .* ma_ratios)
                        elseif outputs[k] == :dudt_grid
                            sum(fic.p_dudt .* ma_ratios)
                        elseif outputs[k] == :dvdt_grid
                            sum(fic.p_dvdt .* ma_ratios)
                        elseif outputs[k] == :si_frac_grid
                            area_tot/LG.area(cell_poly)
                        elseif outputs[k] == :overarea_grid
                            sum(fic.overarea)/length(fic)
                        elseif outputs[k] == :mass_grid
                            mass_tot
                        elseif outputs[k] == :area_grid
                            area_tot
                        elseif outputs[k] == :height_grid
                            sum(fic.height .* ma_ratios)
                        elseif outputs[k] == :stress_xx_grid
                            sum([s[1, 1] for s in fic.stress] .* ma_ratios)
                        elseif outputs[k] == :stress_yx_grid
                            sum([s[1, 2] for s in fic.stress] .* ma_ratios)
                        elseif outputs[k] == :stress_xy_grid
                            sum([s[2, 1] for s in fic.stress] .* ma_ratios)
                        elseif outputs[k] == :stress_yy_grid
                            sum([s[2, 2] for s in fic.stress] .* ma_ratios)
                        elseif outputs[k] == :stress_eig_grid
                            xx = sum([s[1, 1] for s in fic.stress] .* ma_ratios)
                            yx = sum([s[1, 2] for s in fic.stress] .* ma_ratios)
                            xy = sum([s[2, 1] for s in fic.stress] .* ma_ratios)
                            yy = sum([s[2, 2] for s in fic.stress] .* ma_ratios)
                            stress = maximum(eigvals([xx yx; xy yy]))
                            if abs(stress) > 1e8
                                stress = 0.0
                            end
                            stress
                        elseif outputs[k] == :strain_ux_grid
                            sum([s[1, 1] for s in fic.strain] .* ma_ratios)
                        elseif outputs[k] == :strain_vx_grid
                            sum([s[1, 2] for s in fic.strain] .* ma_ratios)
                        elseif outputs[k] == :strain_uy_grid
                            sum([s[2, 1] for s in fic.strain] .* ma_ratios)
                        elseif outputs[k] == :strain_vy_grid
                            sum([s[2, 2] for s in fic.strain] .* ma_ratios)
                        end
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
    getattrs(output::FloeOutput)

Returns unit and comment attributes for each output type to be saved within output NetCDF file
Input:
        output<FloeOutput>
Output:
        <Tuple(String, String)> tuple of string units and comments to be saved to output NetCDF file
"""
function getattrs(output::Symbol)
    unit, comment =
        output == :sim ? ("Struct", "Starting state of simulation after initialization") :
        # Model structs/struct arrays
        output == :ocean ? ("Struct", "Model's ocean struct at given simulation timesteps") :
        output == :atmos ? ("Struct", "Model's atmosphere struct at given simulation timesteps") :
        output == :floes ? ("StructArray", "Model's array of floe structs at given simulation timesteps") :
        # Floe fields
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
        output == :id ? ("unitless", "Floe ID - represents unique floe (ghosts have the same as parents since they are copies)") :
        output == :ghost_id ? ("unitless", "Ghost ID - represents individual ghosts so every ghost of 1 parent has unique ID (parents have ghost ID of 0)") :
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
        output == :stress ? ("N/m^2", "Stress on the floe in the form [xx yx; xy, yy]") :
        output == :stress_history ? ("N/m^2", "List of last stresses felt on the floe from previous timesteps") :
        output == :strain ? ("unitless", "Strain on the floe in the form [ux vx; uy vy]") :
        output == :p_dxdt ? ("m/s", "Floe x-velocity from previous time step") :
        output == :p_dydt ? ("m/s", "Floe y-velocity from previous time step") :
        output == :p_dudt ? ("m/s^2", "Floe x-acceleration from previous time step") :
        output == :p_dvdt ? ("m/s^2", "Floe y-acceleration from previous time step") :
        output == :p_dξdt ? ("rad/s^2", "Floe angular acceleration from the previous timestep") :
        output == :p_dαdt ? ("rad/s", "Floe angular velocity from previous timestep") :
        # Gridded data fields
        output == :u_grid ? ("m/s", "Average floe x-direction velocity in grid cell") :
        output == :v_grid ? ("m/s", "Average floe y-direction velocity in grid cell") :
        output == :dudt_grid ? ("m/s^2", "Average floe x-direction acceleration in grid cell") :
        output == :dvdt_grid ? ("m/s^2", "Average floe y-direction acceleration in grid cell") :
        output == :overarea_grid ? ("m", "Average overlap area of floe with other floes in grid cell") :
        output == :mass_grid ? ("kg", "Average mass of floes in grid cell") :
        output == :area_grid ? ("m^2", "Average area of floes in grid cell") :
        output == :height_grid ? ("m", "Average height of floes in grid cell") :
        output == :si_frac_grid ? ("unitless", "Fraction of grid cell covered by floes") :
        output == :stress_xx_grid ? ("N/m^2", "Average xx stress on floes in a given grid cell") :
        output == :stress_yx_grid ? ("N/m^2", "Average yx stress on floes in a given grid cell") :
        output == :stress_xy_grid ? ("N/m^2", "Average xy stress on floes in a given grid cell") :
        output == :stress_yy_grid ? ("N/m^2", "Average yy stress on floes in a given grid cell") :
        output == :stress_eig_grid ? ("N/m^2", "Maximum eigenvalue of the stress matricies [xx yx; xy yy]") :
        output == :strain_ux_grid ? ("unitless", "Average ux strain on floes in a grid cell") :
        output == :strain_vx_grid ? ("unitless", "Average vx strain on floes in a grid cell") :
        output == :strain_uy_grid ? ("unitless", "Average uy strain on floes in a grid cell") :
        output == :strain_vy_grid ? ("unitless", "Average vy strain on floes in a grid cell") :
        ("", "") # if symbol isn't found, return empty attributes
    return unit, comment
end
