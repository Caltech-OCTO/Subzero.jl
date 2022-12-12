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

"""
    domain_xycoords(domain::Domain)

X-y coordinates for a rectangle in shape of given domain to be used for plotting.
Inputs:
        domain <Domain>
Outputs:
        x-coordinates <Float Vector>
        y-coordinates <Float Vector>
"""
function domain_xycoords(domain::Domain)
    north = domain.north.val
    south = domain.south.val
    east = domain.east.val
    west = domain.west.val
    return [west, west, east, east, west], [north, south, south, north, north]
end

"""
    setup_plot(model::Model)

Plots grid lines, model domain, and model topology on a plot with a ratio determined by the grid coordinates.
This plot will be used throughout simulation so these qualities only need to be plotted/set once as they cannot change during the simulation.

Inputs:
        model <Model>
Output:
        plt <Plots.Plot> with gridlines and aspect ratio set and model domain and topology plotted
"""
function setup_plot(model::Model)
    # Plot Grid Lines and set image ratio
    xmax = model.grid.xg[end]/1000  # kilometers
    ymax = model.grid.yg[end]/1000 # kilometers
    ratio = ymax/xmax
    plt = Plots.plot(xlims = (model.grid.xg[1]/1000, xmax),
                     ylims = (model.grid.yg[1]/1000, ymax),
                     size = (1500, 1200),
                     aspect_ratio=ratio,
                     xlabel = "[km]",
                     ylabel = "[km]")
    # Plot Domain Border
    domainx, domainy = domain_xycoords(model.domain)
    plot!(plt, domainx./1000, domainy./1000 , seriestype = [:shape],
          linecolor = :black, fillalpha = 0.0, lw = 2, legend=false)

    # Plot Topology
    if length(size(model.topos)) > 0
        topo_coords = model.topos.coords
        plot!(plt, [LG.Polygon([c[1] ./ 1000]) for c in topo_coords],
              fill = :grey)
    end
    return plt
end

function setup_plot(domain_data::String, plot_size = (1500, 1200))
    # Open file to get needed values
    plt = jldopen(domain_data, "r") do file
        d = file["domain"]
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
                         ylabel = "[km]")
        if !isempty(d.topography)
            topo_coords = d.topography.coords
            plot!(plt, [LG.Polygon([c[1] ./ 1000]) for c in topo_coords], fillcolor = :grey)
        end
    end
    return plt
end

"""
    plot_sim(model, plt)

Plot ocean velocities and floes for current model time step using given plot.
    
Inputs:
        model   <Model>
        plt     <Plots.plt>
        time    <Int> model timestep
Outpits:
        Save new plot
"""
function plot_sim(model, plt, time)
    # Plot Ocean Vector Field - also clears previous plot
    xgrid, ygrid = Subzero.grids_from_lines(model.grid.xc, model.grid.xc)
    plt_new = quiver(plt, vec(xgrid ./ 1000), vec(ygrid ./ 1000),
            quiver=(vec(model.ocean.u), vec(model.ocean.v)), color = :lightgrey,
            title = string("Time: ", round(time/6, digits = 2), " minutes"))

    # Plot Floes --> only plot "alive" floes
    floe_coords = model.floes.coords
    floe_alive = model.floes.alive
    plot!(plt_new, [LG.Polygon([floe_coords[i][1] ./ 1000]) for i in eachindex(floe_coords) if floe_alive[i] == 1], fill = :lightblue)
          
    # Save plot
    Plots.savefig(plt_new, "figs/collisions/plot_$time.png")
end

function create_sim_gif(floes_fn, domian_fn)
    sim_data = NCDataset(floes_fn)  # TODO: Change this if we don't want NetCDFs

    xcoords = sim_data["xcoords"][:, :, :]
    ycoords = sim_data["ycoords"][:, :, :]
    alive = sim_data["alive"][:, :]
    anim = @animate for tstep in eachindex(sim_data["time"][:])
        plt = setup_plot(domian_fn)  # Uses JLD2
        for i in eachindex(sim_data["floes"][:])
            if alive[tstep, i] == 1
            plot!(plt, xcoords[tstep, i, :]./1000, ycoords[tstep, i, :]./1000,
                           seriestype = [:shape,], fill = :lightblue, legend=false)
            end
        end
    end
    gif(anim, string("figs/collisions/f.gif"), fps = 15)
end