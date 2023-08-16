"""
Plotting functions for Subzero Simulation
"""

"""
    prettytime(t)

Turn time in seconds into units of minutes, hours, days, or years as appropriate
Input:
    t   <Real> number of seconds
Output:
    String with value and units
Note:
    modified from
    https://github.com/JuliaCI/BenchmarkTools.jl/blob/master/src/trials.jl
"""
function prettytime(t)
    minute = 60
    hour = 3600
    day = 24*3600
    year = 360*day

    iszero(t) && return "0 seconds"
    if t < minute
        value = floor(Int, t)
        units = value == 1 ? "second" : "seconds"
    elseif t < hour
        value = floor(Int, t / minute)
        units = value == 1 ? "minute" : "minutes"
    elseif t < day
        value = floor(Int, t / hour)
        units = value == 1 ? "hour" : "hours"
    elseif t < year
        value = floor(Int, t / day)
        units = value == 1 ? "day" : "days"
    else
        value = floor(Int, t / year)
        units = value == 1 ? "year" : "years"
    end
    return @sprintf("%d %s", value, units)
end

"""
    get_curl(fldx,fldy,dx,dy)

Calculate curl from ocean u and v velocity fields
Inputs:
    fldx    <Matrix{AbstractFloat}> ocean u velocity field
    fldy    <Matrix{AbstractFloat}> ocean v velocity field
    dx      <AbstractFloat> x-distance over which velocity fields are provided
    dy      <AbstractFloat> y-distance over which velocity fields are provided
Output:
    <Matrix{AbstractFloat}> curl
"""
function get_curl(fldx,fldy,dx,dy)
    # fldx must be on the u-grid point and fldy on the v-grid
    # returns field on the omega grid
    nx,ny = size(fldx)
    dvdx = zeros(ny,nx)
    dudy = zeros(ny,nx)
    dvdx[1:end-1,:] = diff(fldy,dims=1)/dx
    dudy[:,1:end-1] = diff(fldx,dims=2)/dy
    return (dvdx - dudy)
end

"""
    calc_ro_field(ocean_fn)

Calculate surface vorticity for ocean file.
Inputs:
    ocean_fn    <String> filename for ocean NetCDF filename
Outputs:
    ro  <Array> 3D array where first two dimensions are ocean size (Nx, Ny) and the third
            dimension is time over the simulation
    xc  <Vector> x grid points for ro values
    yc  <Vector> y grid points for ro values
"""
function calc_ro_field(ocean_fn)
    xc = NetCDF.ncread(ocean_fn, "xC")
    yc = NetCDF.ncread(ocean_fn, "yC")
    dx = xc[2] - xc[1]
    dy = yc[2] - yc[1]
    Nx = length(xc)
    Ny = length(yc)
    usurf = NetCDF.ncread(ocean_fn, "u")[1:Nx,1:Ny,:,:]
    vsurf = NetCDF.ncread(ocean_fn, "v")[1:Nx,1:Ny,:,:]
    omega = 2*π/(3600*24)
    f = f = 2*omega*sin(70*π/180)
    nsteps = size(usurf, 4)
    ro = zeros(Nx, Ny, nsteps)
    @views for i in 1:nsteps
        ro[:, :, i] .= get_curl(usurf[:,:,1,i], vsurf[:,:,1,i], dx,dy) ./ f
    end
    return ro, xc, yc
end

"""
    CoordPlot

Recipe for plotting list of PolyVecs that creates two new function coordplot and
coordplot!
"""
@recipe(CoordPlot, coord_list) do scene
    Attributes(
        color = :lightblue,    # floe fill color (could give transparent color)
        strokecolor = :black,  # outline color of floes
        strokewidth = 1,       # width of floe outline
    )
end

"""
    Makie.plot!(coordplot)

Defines coordplot and coordplot! for plotting lists of PolyVecs.
"""
function Makie.plot!(coordplot::CoordPlot{<:Tuple{<:Vector{<:PolyVec}}})
    coord_list = coordplot[1]
    poly_list = @lift([[Point2f(verts) for verts in c[1]] for c in $coord_list])
    Makie.poly!(
        coordplot,
        poly_list,
        color = coordplot[:color],
        strokecolor = coordplot[:strokecolor],
        strokewidth = coordplot[:strokewidth],    
    )
    coordplot
end

"""
    plot_sim(
        floe_fn,
        initial_state_fn,
        title,
        Δt,
        output_fn
    )

Basic plotting of sim as an example of how to create video. Does not have
underlying ocean.
Inputs:
    floe_fn           <String> Subzero floe outputwriter output file
    initial_state_fn  <String> Subzero initial state writer output file
    title             <String> plot title
    Δt                <Int> length of timestep in seconds
    output_fn         <String> output filename (should be .mp4)
Output:
    Saves video as output_fn
"""
function plot_sim(
    floe_fn,
    initial_state_fn,
    Δt,
    output_fn,
)
    # Open files
    file = jldopen(floe_fn)
    domain = load(initial_state_fn)["sim"].model.domain
    timesteps = keys(file["centroid"])
    # Set up observables
    floes = Observable(file["coords"][timesteps[1]])
    # Plot floes
    fig, ax, _ = coordplot(floes)
    # Set axis limits and names
    xlims!(domain.west.val, domain.east.val)
    ylims!(domain.south.val, domain.north.val)
    ax.xlabel =  "Meters"
    ax.ylabel = "Meters"
    # Plot topography
    if !isempty(domain.topography)
        coordplot!(domain.topography.coords, color = :lightgrey)
    end
    # Create movie
    record(fig, output_fn, timesteps; framerate = 20) do time
        ax.title = Subzero.prettytime(parse(Float64, time) * Δt)
        new_coords = file["coords"][time]
        floes[] = new_coords
    end
    close(file)
end


"""
    plot_sim_with_ocean_field(
        floe_fn,
        initial_state_fn,
        Δt,
        ocean_fn,
        ocean_func,
        output_fn,
    )

Basic plotting of sim as an example of how to create video. Does not have
underlying ocean.
Inputs:
    floe_fn           <String> Subzero floe outputwriter output file
    initial_state_fn  <String> Subzero initial state writer output file
    Δt                <Int> length of timestep in seconds
    ocean_fn          <String> Oceananigans output surface.nc file
    ocean_func        <Function> function that takes in ocean_fn and returns a
                        Nx by Ny by timesteps surface field for plotting as well
                        as xc and yc fields (see calc_ro_field for example)
    colorbar_title    <String> name for colorbar associated with ocean_func vals
    output_fn         <String> output filename (should be .mp4)
Output:
    Saves video as output_fn.
"""
function plot_sim_with_ocean_field(
    floe_fn,
    initial_state_fn,
    Δt,
    ocean_fn,
    ocean_func,
    colorbar_title,
    output_fn,
)
    # Open files
    file = jldopen(floe_fn)
    domain = load(initial_state_fn)["sim"].model.domain
    timesteps = keys(file["centroid"])
    ocean_data, xc, yc = ocean_func(ocean_fn)
    # Set up observables needed for plotting
    floes = Observable(file["coords"][timesteps[1]])
    ocean_vals = Observable(@view ocean_data[:, :, 1])
    min_ocn_val, max_ocn_val = extrema(ocean_data)
    fig = Figure()
    # Plot ocean
    ax, hm = heatmap(
        fig[1, 1],
        xc,
        yc,
        ocean_vals,
        colormap = :RdBu_9,
        colorrange = (min_ocn_val, max_ocn_val)
    )
    # Add axis limits and titles
    xlims!(domain.west.val, domain.east.val)
    ylims!(domain.south.val, domain.north.val)
    ax.xlabel =  "Meters"
    ax.ylabel = "Meters"
    # Add colorbar
    Colorbar(fig[1, 2], hm, label = colorbar_title)
    # Plot floes
    coordplot!(fig[1, 1], floes)
    # Plot topography
    if !isempty(domain.topography)
        coordplot!(fig[1, 1], domain.topography.coords, color = :lightgrey)
    end

    # Create movie
    record(fig, output_fn, 1:length(timesteps), framerate = 20) do i
        time = timesteps[i]
        ax.title = prettytime(parse(Float64, time) * Δt)
        new_coords = file["coords"][time]
        ocean_vals[] = @view ocean_data[:, :, i]
        floes[] = new_coords
    end
    close(file)
end