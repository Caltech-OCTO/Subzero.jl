module SubzeroMakieExt

using CairoMakie
using Subzero, JLD2
import Subzero: prettytime, plot_sim, plot_sim_with_ocean_field

"""
    CoordPlot

Recipe for plotting list of `PolyVec`s that creates two new functions: `coordplot` and
`coordplot!`.

These functions each take in a vector of `PolyVec`s to plot the floe field at a given
timestep.

Note: I would like to remove this recipe and it's corresponding `plot!` function eventually.
Now that all of the floes' carry around a `GeoInterface` polygon, we can simply use
`GeoInterfaceMakie` to plot them and not deal with plotting from coordiantes.
"""
CairoMakie.@recipe(CoordPlot, coord_list) do scene
    Attributes(
        color = :lightblue,    # floe fill color (could give transparent color)
        strokecolor = :black,  # outline color of floes
        strokewidth = 1,       # width of floe outline
    )
end

"""
    Makie.plot!(coordplot)

Defines coordplot and coordplot! for plotting vectors of PolyVecs, representing the floe
field at a given timestep.

Note: I would like to remove this function and it's corresponding `CoordPlot` recipe
eventually. Now that all of the floes' carry around a `GeoInterface` polygon, we can simply
use `GeoInterfaceMakie` to plot them and not deal with plotting from coordiantes.
"""
function CairoMakie.plot!(coordplot::CoordPlot{<:Tuple{<:Vector{<:PolyVec}}})
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
    plot_sim(floe_fn, initial_state_fn, title, Δt, output_fn)

Basic plotting of a simulation using the simulation's floe and initial state files. This
function is meant for basic plotting and as an example of how to create video.
Does not have underlying ocean.

## Arguments:
- `floe_fn::String`: $(Subzero.FLOE_FN_DEF)
- `initial_state_fn::String`: $(Subzero.INITIAL_STATE_FN_DEF)
- `title::String`: plot title
- `Δt::Int`: $(Subzero.ΔT_DEF)
- `output_fn::String`: $(Subzero.MP4_OUTPUT_FN)
"""
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

"""
    plot_sim_with_ocean_field(floe_fn, initial_state_fn, Δt, ocean_fn, ocean_func, output_fn)

Basic plotting of sim as an example of how to create video. Does not have
underlying ocean.

Basic plotting of a simulation using the simulation's floe and initial state files, as well
as an surface data file output by `Oceananigans` that can be run through an `ocean_func` to
get a gridded output. This function is meant for basic plotting and as an example of how to
create video.

## Arguments:
- `floe_fn::String`: $(Subzero.FLOE_FN_DEF)
- `initial_state_fn::String`: $(Subzero.INITIAL_STATE_FN_DEF)
- `Δt::Int`: $(Subzero.ΔT_DEF)
- `ocean_fn::String`: $(Subzero.OCEAN_FN_DEF)
- `ocean_func::Function`: function that takes in `ocean_fn` and returns an Nx by Ny by timesteps surface field for plotting as well as xc and yc fields (see `calc_ro_field` for example)
- `colorbar_title::String`: name for colorbar associated with ocean_func vals
- `output_fn::String`: $(Subzero.MP4_OUTPUT_FN)
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

end  # module