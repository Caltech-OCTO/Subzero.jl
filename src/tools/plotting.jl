"""
Plotting functions for Subzero Simulation
"""

"""
grids_from_lines(xlines, ylines)

Creates x-grid and y-grid. Assume xlines has length n and ylines has length m.
xgrid is the grid's xline vector repeated m times as rows in a mxn array and
ygrid is the yline vector repeated n times as columns in a mxn vector. xlines
and ylines are typically either xg and yg or xc and yc.
"""
function grids_from_lines(xlines, ylines)
    xgrid = repeat(reshape(xlines, 1, :), inner=(length(ylines),1))
    ygrid = repeat(ylines, outer = (1, length(xlines)))
    return xgrid, ygrid
end

#------------------ Plotting from FloeWriter Files Post-Simulation ------------------#

"""
    setup_plot(domain_data::String, plot_size = (1500, 1200))

Set up plot object using domain data from file. x and y limits based on domain
and topograhy plotted in grey.
Inputs:
    domain_fn  <String> file path to JLD2 file holding domain struct information
    plot_size  <Tuple{Int, Int}> size of output plot - default is (1500, 1200)
    plot_ocn   <Boolean> boolean flag to plot the ocean arrows if true
Outputs:
    Plot with x and y xlim determed by domain and including all topography. 
"""
function setup_plot(init_pn::String, plot_size = (1500, 1500), plot_ocn = false)
    # Open file to get needed values
    file = jldopen(init_pn, "r")
    g = file["sim"].model.grid
    d = file["sim"].model.domain
    o = file["sim"].model.ocean
    
    xmin = d.west.val/1000
    xmax = d.east.val/1000
    ymin = d.south.val/1000
    ymax = d.north.val/1000

    # Plot domain using file data
    ratio = (ymax-ymin)/(xmax-xmin)
    plt = Plots.plot(xlims = (xmin, xmax),
                        ylims = (ymin, ymax),
                        size = plot_size,
                        aspect_ratio=ratio,
                        xlabel = "[km]",
                        ylabel = "[km]",
                        xtickfontsize = 20,
                        ytickfontsize = 20,
                        xguidefontsize=25,
                        yguidefontsize=25,
                        margin = 10mm,
                        )
    if !isempty(d.topography)
        topo_coords = separate_xy.(d.topography.coords)
        plot!(
            plt,
            first.(topo_coords)./1000,
            last.(topo_coords)./1000,
            seriestype = [:shape,],
            fill = :grey,
            legend=false,
        )
    end
    if plot_ocn
        xg = g.x0:g.Δx:g.xf
        yg = g.y0:g.Δy:g.yf
        xgrid, ygrid = Subzero.grids_from_lines(xg, yg)
        quiver!(
            plt,
            vec(xgrid ./ 1000),
            vec(ygrid ./ 1000),
            quiver=(
                vec(o.u'),
                vec(o.v'),
            ),
            color = :lightgrey,
        )
    end
    JLD2.close(file)
    return plt
end

"""
    create_sim_gif(floes_pn, domain_fn, output_fn, plot_size = (1500, 1500))

Create a gif of a simulation given a file with domain information and floe
information from a floe output writer.
Inputs:
    floes_pn   <String> file path to file output by floe output writer that
                    inclues coordinate and status fields
    init_pn    <String> file path to JLD2 file holding domain struct information
    output_fn  <String> file path to save gif
    plot_size  <Tuple(Int, Int)> size of output gif in pixels - default
                    (1500, 1500)
    fps        <Int> frames per second
    plot_ocn   <Boolean> boolean flag to plot the ocean arrows if true
Outputs:
    Saves simulation gif with floes and topography plotted.
"""
function create_sim_gif(
    floe_pn,
    init_pn,
    output_fn;
    plot_size = (1500, 1500),
    fps = 15,
    plot_ocn = false,
)
    # Get floe data
    floe_data = jldopen(floe_pn, "r")
    status = floe_data["status"]
    coords = floe_data["coords"]
    times = keys(coords)
    # Plot floe data
    anim = @animate for t in times
        plt = setup_plot(init_pn, plot_size, plot_ocn)
        verts = Subzero.separate_xy.(coords[t])
        for i in eachindex(verts)
            if status[t][i].tag != remove
                plot!(
                    plt,
                    first(verts[i])./1000,
                    last(verts[i])./1000,
                    seriestype = [:shape,],
                    fill = :lightblue,
                    legend=false,
                )
            end
        end
    end
    JLD2.close(floe_data)
    gif(anim, output_fn, fps = fps)
    return
end

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
    dvdx = zeros((ny,nx))
    dudy = zeros((ny,nx))
    dvdx[1:end-1,:] = diff(fldy,dims=1)/dx
    dudy[:,1:end-1] = diff(fldx,dims=2)/dy
    return (dvdx - dudy)
end

"""
    create_coupled_ro_sim_gif(
        floe_pn,
        init_pn,
        surface_pn,
        output_fn;
        plot_size = (1500, 1500),
        fps = 15,
    )

Plots floes over Oceananigans surface Ro = ξ/f.
Inputs:
    floes_pn    <String> file path to file output by floe output writer that
                    inclues coordinate and status fields
    init_pn     <String> file path to JLD2 file holding domain information
    surface_pn  <String> file path to NetCDF file holding Oceananigans surface
                    data
    output_fn   <String> file path to save gif
    plot_size   <Tuple(Int, Int)> size of output gif in pixels - default
                    (1500, 1500)
    fps         <Int> frames per second - default 15
Output:
    Saves gif to output_fn
"""
function create_coupled_ro_sim_gif(
    floe_pn,
    init_pn,
    surface_pn,
    output_fn;
    plot_size = (1500, 1500),
    fps = 15,
)
    # Loading surface data
    fname = surface_pn
    xc = NetCDF.ncread(fname, "xC")
    yc = NetCDF.ncread(fname, "yC")
    Nx = length(xc)
    Ny = length(yc)
    usurf = NetCDF.ncread(fname, "u")[1:Nx,1:Ny,:,:]
    vsurf = NetCDF.ncread(fname, "v")[1:Nx,1:Ny,:,:]
    t = NetCDF.ncread(fname, "time")
    nit = length(t)
    dx = xc[2] - xc[1]
    dy = yc[2] - yc[1]
    omega = 2*π/(3600*24)
    f = 2*omega*sin(70*π/180)

    # Get floe data
    floe_data = jldopen(floe_pn, "r")
    status = floe_data["status"]
    coords = floe_data["coords"]

    # Plot floe data
    plt = setup_plot(init_pn, plot_size)
    times = keys(coords)
    anim = @animate for i = 1:nit
        ts = times[i]
        # get floe coords at timestep
        verts = Subzero.separate_xy.(coords[ts])
        # calculate Ro for timestep
        vort = get_curl(usurf[:,:,1,i],vsurf[:,:,1,i],dx,dy) # curl
        Ro = vort/f
        # Setup plot for timestep
        plt = setup_plot(init_pn, plot_size)
        # Plot ocean for timestep
        contour!(
            plt,
            xc/1000,
            yc/1000,
            Ro',  # should be transposed 
            xlabel="x [km]",
            ylabel="y [km]",
            title=prettytime(t[i]),
            titlefontsize = 35,
            linewidth = 0,
            fill=true,
            color=:balance,
            clims=(-0.6, 0.6),
            legend=false,
            colorbar=true,
            colorbar_title = "Ro = ζ/f",
            colorbar_titlefontsize = 25,
        )
        # Plot floes for timestep
        for j in eachindex(verts)
            if status[ts][j].tag != remove
                plot!(
                    plt,
                    first(verts[j])./1000,
                    last(verts[j])./1000,
                    seriestype = [:shape,],
                    fill = :lightblue,
                    legend=false,
                )
            end
        end
    end
    JLD2.close(floe_data)
    gif(anim, output_fn, fps = fps)
    return
end

# using Makie, CairoMakie
# Makie.@recipe(FloeScene) do scene
#     Attributes(
#         floecolor = :lightblue,
#         oceancolor = :gray,
#     )
# end

# function Makie.plot!(
#     fs::FloeScene{
#         <:Tuple{
#             <:Real # timestep
#             <:AbstractMatrix{<:Real},  # floe coordinates 2xn
#             <:AbstractMatrix{<:Real},  # ocean x
#             <:AbstractMatrix{<:Real},  # ocean y
#             <:AbstractMatrix{<:Real},  # ocean x tracer/vector
#             <:AbstractMatrix{<:Real},  # ocean y tracer/vector
#         },
#     },
# )
#     timestep = fs[1]
#     floecoords = fs[2]

#     points = Observable(Point2f[])
#     function update_plot(timestep, floecoords)
#         empty!(points[])
#         for i in eachrow(floecoords)
#             push!(points, Point2f(floecoords[i, 1], floecoords[i, 2]))
#         end
#     end
#     Makie.Observables.onany(update_plot, timestep, floecoords)
#     update_plot(timestep[], floecoords[])

#     title!(fs, string(timestep))

#     poly!(fs, points, color = fs.floecolor)

#     return fs
# end

#------------ Plotting for Debugging During Simulation Run --------------#

"""
    setup_plot(model::Model)

Plots grid lines, model domain, and model topology on a plot with a ratio
determined by the grid coordinates. Used for plotting during simulation, mainly
for debugging purposes. These qualities only need to be plotted/set once as
they cannot change during the simulation.

Inputs:
    model <Model>
Output:
    plt <Plots.Plot> with gridlines and aspect ratio set and model domain and
            topology plotted
"""
function setup_plot(model::Model)
    # Plot Grid Lines and set image ratio
    xmin = model.domain.west.val/1000
    xmax = model.domain.east.val/1000
    ymin = model.domain.south.val/1000
    ymax = model.domain.north.val/1000

    # Plot domain using file data
    ratio = (ymax-ymin)/(xmax-xmin)
    plt = Plots.plot(xlims = (xmin, xmax),
                        ylims = (ymin, ymax),
                        size = plot_size,
                        aspect_ratio=ratio,
                        xlabel = "[km]",
                        ylabel = "[km]")

    # Plot Topology
    if !isempty(model.domain.topography)
        topo_coords = model.domain.topography.coords
        plot!(
            plt,
            [LG.Polygon([c[1] ./ 1000]) for c in topo_coords],
            fill = :grey,
        )
    end
    return plt
end

"""
    plot_sim(model, plt)

Plot ocean velocities and floes for current model time step using given plot.
Used for plotting during simulation, mainly for debugging purposes. 
    
Inputs:
        model   <Model>
        plt     <Plots.plt>
        time    <Int> model timestep
Outputs:
        Save new plot.
"""
function plot_sim_timestep(model, plt, time)
    # Plot Ocean Vector Field - also clears previous plot
    xgrid, ygrid = Subzero.grids_from_lines(model.grid.xc, model.grid.xc)
    plt_new = quiver(plt, vec(xgrid ./ 1000), vec(ygrid ./ 1000),
            quiver=(vec(model.ocean.u), vec(model.ocean.v)), color = :lightgrey,
            title = string("Time: ", round(time/6, digits = 2), " minutes"))

    # Plot Floes --> only plot "active" floes
    floe_coords = model.floes.coords
    floe_status = model.floes.status
    plot!(
        plt_new,
        [LG.Polygon(
            [floe_coords[i][1] ./ 1000]
        ) for i in eachindex(floe_coords) if floe_status[i].tag != remove],
        fill = :lightblue,
    )
          
    # Save plot
    Plots.savefig(plt_new, "figs/collisions/plot_$time.png")
end