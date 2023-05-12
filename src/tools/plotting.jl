"""
Plotting functions for Subzero Simulation
"""

"""
grids_from_lines(xlines, ylines)

Creates x-grid and y-grid. Assume xlines has length n and ylines has length m. xgrid is the grid's xline vector repeated m times as rows in a mxn array and ygrid is the yline vector repeated n times as columns in a mxn vector. xlines and ylines are typically either xg and yg of xc and yc.
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
Outputs:
    Plot with x and y xlim determed by domain and including all topography. 
"""
function setup_plot(init_pn::String, plot_size = (1500, 1500))
    # Open file to get needed values
    file = jldopen(init_pn, "r")
    d = file["sim"].model.domain
    
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
Outputs:
    Saves simulation gif with floes and topography plotted.
"""
function create_sim_gif(
    floe_pn,
    init_pn,
    output_fn;
    plot_size = (1500, 1500),
    fps = 15,
)
    # Get floe data
    floe_data = jldopen(floe_pn, "r")
    status = floe_data["status"]
    coords = floe_data["coords"]
    # Plot floe data
    plt = setup_plot(init_pn, plot_size)
    times = keys(coords)
    anim = @animate for t in times
        new_frame = plot(plt)
        verts = Subzero.separate_xy.(coords[t])
        for i in eachindex(verts)
            if status[t][i].tag != remove
                plot!(
                    new_frame,
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

# function prettytime(t)
#     # Modified from: https://github.com/JuliaCI/BenchmarkTools.jl/blob/master/src/trials.jl
#     iszero(t) && return "0 seconds"

#     t = maybe_int(t)

#     if t < 1e-9
#         # No point going to picoseconds, just return something readable.
#         return @sprintf("%.3e seconds", t)
#     elseif t < 1e-6
#         value, units = t * 1e9, "ns"
#     elseif t < 1e-3
#         value, units = t * 1e6, "μs"
#     elseif t < 1
#         value, units = t * 1e3, "ms"
#     elseif t < minute
#         value = t
#         units = value == 1 ? "second" : "seconds"
#     elseif t < hour
#         value = maybe_int(t / minute)
#         units = value == 1 ? "minute" : "minutes"
#     elseif t < day
#         value = maybe_int(t / hour)
#         units = value == 1 ? "hour" : "hours"
#     elseif t < year
#         value = maybe_int(t / day)
#         units = value == 1 ? "day" : "days"
#     else
#         value = maybe_int(t / year)
#         units = value == 1 ? "year" : "years"
#     end

#     if isinteger(value)
#         return @sprintf("%d %s", value, units)
#     else
#         return @sprintf("%.3f %s", value, units)
#     end
# end

# function create_coupled_sim_gif(
#     floe_pn,
#     init_pn,
#     output_fn;
#     plot_size = (1500, 1500),
#     fps = 15,
# )

#     # Loading surface data
#     fname = "../data/" * run_name * "/surface.nc"
#     xc = ncread(fname, "xC")
#     yc = ncread(fname, "yC")
#     Nx = length(xc)
#     Ny = length(yc)
#     usurf = ncread(fname, "u")[1:Nx,1:Ny,:,:]
#     vsurf = ncread(fname, "v")[1:Nx,1:Ny,:,:]
#     wsurf = ncread(fname, "w")
#     Tsurf = ncread(fname, "T")
#     Ssurf = ncread(fname, "S")
#     t = ncread(fname, "time")
#     sigma = gsw_sigma0.(Ssurf, Tsurf)
#     nit = length(t)
#     dx = xc[2] - xc[1]
#     dy = yc[2] - yc[1]
#     Lx = xc[end] + dx/2
#     Ly = yc[end] + dy/2
#     omega = 2*π/(3600*24)
#     f = 2*omega*sin(70*π/180)

#     # Get floe data
#     floe_data = jldopen(floe_pn, "r")
#     status = floe_data["status"]
#     coords = floe_data["coords"]

#     # Plot floe data
#     plt = setup_plot(init_pn, plot_size)
#     times = keys(coords)
#     anim = @animate for t in times
#         new_frame = plot(plt)
#         # plot ocean
#         vort = get_curl(usurf[:,:,1,i],vsurf[:,:,1,i],dx,dy) # curl
#         verts = Subzero.separate_xy.(coords[t])
#         Ro = vort/f
#         conv = -get_div(usurf[:,:,1,i],vsurf[:,:,1,i],dx,dy)
#         title_plt = plot(
#             title = prettytime(t[i]),
#             grid = false,
#             showaxis = false,
#             bottom_margin = -30Plots.px,
#             ticks=nothing,
#         )
#         # plot floes
#         for i in eachindex(verts)
#             if status[t][i].tag != remove
#                 plot!(
#                     new_frame,
#                     first(verts[i])./1000,
#                     last(verts[i])./1000,
#                     seriestype = [:shape,],
#                     fill = :lightblue,
#                     legend=false,
#                 )
#             end
#         end
#     end
#     JLD2.close(floe_data)
#     gif(anim, output_fn, fps = fps)
#     return
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