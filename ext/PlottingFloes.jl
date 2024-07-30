module PlottingFloesExt

using Subzero, Makie, CairoMakie, Printf, GeometryBasics

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

Makie.@recipe(CoordPlot, coord_list) do scene
    Attributes(
        color = :lightblue,    # floe fill color (could give transparent color)
        strokecolor = :black,  # outline color of floes
        strokewidth = 1,       # width of floe outline
    )
end

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

function plot_sim(
    floe_fn,
    initial_state_fn,
    Δt,
    output_fn;
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


end  # module