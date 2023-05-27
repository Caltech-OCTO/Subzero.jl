"""
Structs and functions for fracturing floes
"""

"""
    AbstractFractureCriteria

Abstract type for fracture criteria. Each struct of this type must have a
vertices field representing the criteria in principal stress space. For a given
polygon, the minimum and maximum eigenvalues of its stress field will be its
location in principal stress space. If that stress point falls outside of
the criteria-verticies defined polygon it is a stress great enough to fracture
the floe. Otherwise the floe will not be fractured.
Each fracture criteria type must also have an update_criteria! function defined
that is used to update the criteria each timestep. If the criteria does not need
to be updated, this function can be empty.
"""
abstract type AbstractFractureCriteria end

"""
    NoFracture<:AbstractFractureCriteria

Simple AbstractFractureCriteria type representing when fracturing functionality
is turned off. If this is the type provided to the simulation's FractureSettings
then fractures will not occur.
"""
struct NoFracture<:AbstractFractureCriteria end


"""
    struct HiblerYieldCurve{FT<:AbstractFloat}<:AbstractFractureCriteria
        pstar::FT
        c::FT
        vertices::PolyVec{FT}
    end

Type of AbstractFractureCriteria using TODO: finish! 

Fields:
    pstar       <AbstractFloat> used to tune ellipse for optimal fracturing
    c           <AbstractFloat> used to tune ellipse for optimal fracturing
    verticies   <PolyVec> vertices of criteria in principal stress space
Note:
    Hibler's paper says that: Both pstar and c relate the ice strength to the
    ice thickness and compactness. c is determined to that 10% open water
    reduces the strength substantially and pstar is considered a free parameter
"""
mutable struct HiblerYieldCurve{FT<:AbstractFloat}<:AbstractFractureCriteria
    pstar::FT
    c::FT
    vertices::PolyVec{FT}

    function HiblerYieldCurve{FT}(
        pstar::Real,
        c::Real,
        vertices::PolyVec
    ) where {FT<:AbstractFloat}
        try
            valid_polyvec!(vertices)
        catch
            throw(ArgumentError("The given vertices for the HiblerYieldCurve \
                can't be made into a valid polygon and thus the initial yield \
                curve can't be created."))
        end
        new{FT}(pstar, c, vertices)
    end
end

"""
    calculate_hibler(floes, pstar, c)

Calculate Hibler's Elliptical Yield Curve as described in his 1979 paper
"A Dynamic Thermodynamic Sea Ice Model".
Inputs:
    floes   <StructArray{Floes}> model's list of floes
    pstar   <AbstractFloat> used to tune ellipse for optimal fracturing
    c       <AbstractFloat> used to tune ellipse for optimal fracturing
Outputs:
    vertices <PolyVec{AbstractFloat}> vertices of elliptical yield curve
Note:
    Hibler's paper says that: Both pstar and c relate the ice strength to the
    ice thickness and compactness. c is determined to that 10% open water
    reduces the strength substantially and pstar is considered a free parameter. 
"""
function calculate_hibler(mean_height, pstar, c)
    compactness = 1  # Could be a user input with future development
    p = pstar*mean_height*exp(-c*(1-compactness))
    t = range(0, 2π, length = 100)
    a = p*sqrt(2)/2
    b = a/2
    x = a*cos.(t)
    y = b*sin.(t)
    vertices = [splitdims([x'; y'])]
    rotate_degrees!(vertices, 45)
    translate!(vertices, -p/2, -p/2)
    return valid_polyvec!(vertices)
end

"""
    HiblerYieldCurve(floes, pstar = 2.25e5, c = 20.0)

Calculates Hibler's Elliptical Yield curve using parameters pstar, c, and the
current floe field. 
Inputs:
    floes   <StructArray{Floes}> model's list of floes
    pstar   <AbstractFloat> used to tune ellipse for optimal fracturing
    c       <AbstractFloat> used to tune ellipse for optimal fracturing
Outputs:
    HiblerYieldCurve struct with vertices determined using the calculate_hibler
    function.
"""
HiblerYieldCurve{FT}(
    floes,
    pstar = 2.25e5,
    c = 20.0,
) where {FT <: AbstractFloat}=
    HiblerYieldCurve{FT}(
        pstar,
        c,
        calculate_hibler(mean(floes.height), pstar, c),
    )

"""
    HiblerYieldCurve(args...; kwargs...)

If a float type isn't specified, HiblerYieldCurve will be Float64. Use 
HiblerYieldCurve{Float32}(args...) for HiblerYieldCurve with type
Float32.
"""
HiblerYieldCurve(args...) = HiblerYieldCurve{Float64}(args...)

"""
    update_criteria!(criteria::HiblerYieldCurve, floes)

Update the Hibler Yield Curve vertices based on the current set of floes. The
criteria changes based off of the average height of the floes.
Inputs:
    criteria    <HiblerYieldCurve> simulation's fracture criteria
    floes       <StructArray{Floe}> model's list of floes
Outputs:
    None. Updates the criteria's vertices field to update new criteria. 
"""
function update_criteria!(criteria::HiblerYieldCurve, floes)
    criteria.vertices = calculate_hibler(
        mean(floes.height),
        criteria.pstar,
        criteria.c
    )
end

"""
    determine_fractures(
        floes,
        criteria,
        min_floe_area,
    )

Determines which floes will fracture depending on the principal stress criteria.
Inputs:
    floes           <StructArray{Floe}> model's list of floes
    criteria        <AbstractFractureCriteria> fracture criteria
    min_floe_area   <AbstractFloat> minimum floe area - floes under this area
                        will not be fractured
Outputs:
    <Vector{Int}> list of indices of floes to fracture 
"""
function determine_fractures(
    floes,
    criteria::AbstractFractureCriteria,
    min_floe_area,
)
    # Determine if floe stresses are in or out of criteria allowable regions
    update_criteria!(criteria, floes)
    σvals = combinedims(sort.(eigvals.(floes.stress)))'
    in_idx = points_in_poly(σvals, criteria.vertices)
    # If stresses are outside of criteria regions, we will fracture the floe
    frac_idx = .!(in_idx)
    frac_idx[floes.area .< min_floe_area] .= false
    return range(1, length(floes))[frac_idx]
end

"""
    deform_floe!(
        floe,
        deformer_coords,
        deforming_forces,
    )

Deform a floe around the area of its collision with largest area overlap within
the last timestep.
Inputs:
        floe                <Floe> floe to deform
        deformer_coords     <PolyVec> coords of floe that is deforming floe
                                argument
        deforming_forces    <Vector{AbstractFloat}> 1x2 matrix of forces between
                                floe and the deforming floe from floe's
                                interactions - of the form: [xforce yforce] 
Outputs:
        None. The input floe's centroid, coordinates, and area are updated to
        reflect a deformation due to the collision with the deforming floe. 
"""
function deform_floe!(
    floe,
    deformer_coords,
    deforming_forces,
    consts,
    Δt,
    mc_n,
    rng,
)
    poly = LG.Polygon(floe.coords)
    deformer_poly = LG.Polygon(deformer_coords)
    overlap_region = sortregions(LG.intersection(poly, deformer_poly))[1]
    # If floe and the deformer floe have an overlap area
    if LG.area(overlap_region) > 0
        # Determine displacement of deformer floe
        rcent = find_poly_centroid(overlap_region)
        dist = calc_point_poly_dist(
            [rcent[1]],
            [rcent[2]],
            find_poly_coords(overlap_region),
        )
        force_fracs = deforming_forces ./ 2norm(deforming_forces)
        Δx, Δy = abs.(dist)[1] .* force_fracs
        # Temporarily move deformer floe to find new shape of floe
        deformer_poly = LG.Polygon(translate(deformer_coords, Δx, Δy))
        new_floe_poly = sortregions(LG.difference(poly, deformer_poly))[1]
        new_floe_area = LG.area(new_floe_poly)
        # If didn't change floe area by more than 90%
        if new_floe_area > 0 && new_floe_area/floe.area > 0.9
            # Update floe shape and conserve mass
            moment_tmp = floe.moment
            x_tmp, y_tmp = floe.centroid
            replace_floe!(
                floe,
                new_floe_poly,
                floe.mass,
                consts,
                mc_n,
                rng,
            )
            conserve_momentum_combination!(
                floe.mass,
                moment_tmp,
                x_tmp,
                y_tmp,
                Δt,
                floe,
            )
        end
    end
    return
end

"""
    split_floe(
        floe,
        npieces,
        rng,
        ::Type{T} = Float64
    )
Splits a given floe into pieces using voronoi tesselation.
User will recieve a warning if floe isn't split.
Inputs:
    floe    <Floe> floe in simulation
    npieces <Int> number of pieces to try to split the floe into - voronoi
                tesselation has an element of randomness so this number is not
                guarenteed but user will be warned if floe isn't split at all
    rng     <RNG> random number generator used for voronoi tesselation
            <Type{T}> AbstractFloat type that used for simulation calculations
Outputs:
    new_floes   <StructArray{Floes}> list of pieces floe is split into, each of
                    which is a new floe
"""
function split_floe(
    floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    rng,
    fracture_settings,
    coupling_settings,
    consts,
    Δt,
) where {FT}
    new_floes = StructArray{Floe{FT}}(undef, 0)
    # Generate voronoi tesselation in floe's bounding box
    scale_fac = fill(2floe.rmax, 2)
    trans_vec = [floe.centroid[1] - floe.rmax, floe.centroid[2] - floe.rmax]
    pieces = generate_voronoi_coords(
        fracture_settings.npieces,
        scale_fac,
        trans_vec,
        [floe.coords],
        rng,
        1,  # Warn if only 1 point is identified as the floe won't be split
    )
    if !isempty(pieces)
        # Intersect voronoi tesselation pieces with floe
        floe_poly = LG.Polygon(floe.coords)
        pieces_polys = [rmholes(
            LG.intersection(LG.Polygon(p), floe_poly)
        ) for p in pieces]
        # Conserve mass within pieces
        pieces_areas = [LG.area(p) for p in pieces_polys]
        total_area = sum(pieces_areas)
        # Create floes out of each piece
        for i in eachindex(pieces_polys)
            if pieces_areas[i] > 0
                mass = floe.mass * (pieces_areas[i]/total_area)
                height = mass / (consts.ρi * pieces_areas[i])
                pieces_floes = poly_to_floes(
                    FT,
                    pieces_polys[i],
                    height,
                    0;  # Δh - range of random height difference between floes
                    ρi = consts.ρi,
                    u = floe.u,
                    v = floe.v,
                    ξ = floe.ξ,
                    mc_n = coupling_settings.mc_n,
                    nhistory = fracture_settings.nhistory,
                    rng = rng,
                )
                append!(new_floes, pieces_floes)
            end
        end
    end

    # Conserve momentum and update strain
    conserve_momentum_fracture!(
        floe,
        new_floes,
        Δt,
    )
    for s in new_floes.strain
        s .= floe.strain
    end

    return new_floes
end

"""
    fracture_floes!(
        floes,
        max_floe_id,
        rng,
        fracture_settings,
        coupling_settings,
        simp_settings,
        consts::Constants{FT},
        Δt,
    )
Fractures floes that meet the criteria defined in the fracture settings.
Inputs:
    floes       <StructArray{Floe}> model's list of floes
    max_floe_id <Int> maximum ID of any floe created so far in simulation
    rng         <RNG> random number generator
    fracture_settings   <FractureSettings> sim's fracture settings
    coupling_settings   <CouplingSettings> sim's coupling settings
    simp_settings       <SimplificationSettings> sim's simplification settings
    consts              <Constants> sim's constants
    Δtout               <Int> length of simulation timestep in seconds
Outputs:
    max_floe_id <Int> new highest floe ID after adding new floes to floe array.
    Floe pieces added to floe array and original fractured floes removed.
"""
function fracture_floes!(
    floes::StructArray{Floe{FT}},
    max_floe_id,
    rng,
    fracture_settings,
    coupling_settings,
    simp_settings,
    consts,
    Δt,
) where {FT <: AbstractFloat}
    # Determine which floes will fracture
    frac_idx = determine_fractures(
        floes,
        fracture_settings.criteria,
        simp_settings.min_floe_area,
    )
    # Initialize list for new floes created from fracturing existing floes
    nfloes2frac = length(frac_idx)
    fracture_list = [StructArray{Floe{FT}}(undef, 0) for _ in 1:nfloes2frac]
    # Fracture floes that meet criteria 
    Threads.@threads for i in 1:nfloes2frac
        # Deform floe around largest impact site
        if fracture_settings.deform_on
            inters = floes.interactions[frac_idx[i]]
            inters = inters[.!(isinf.(inters[:, floeidx])), :]
            if !isempty(inters)
                _, max_inters_idx = findmax(inters[:, overlap])
                deforming_inter = inters[max_inters_idx, :]
                deforming_floe_idx = Int(deforming_inter[floeidx])
                if deforming_floe_idx <= length(floes)
                    deform_floe!(
                        LazyRow(floes, frac_idx[i]), 
                        floes.coords[deforming_floe_idx],
                        deforming_inter[xforce:yforce],
                        consts,
                        Δt,
                        coupling_settings.mc_n,
                        rng,
                    )
                end
            end
        end
        # Split flie into pieces
        new_floes = split_floe(
            LazyRow(floes, frac_idx[i]),
            rng,
            fracture_settings,
            coupling_settings,
            consts,
            Δt,
        )
        append!(fracture_list[i], new_floes)
    end
    # Remove old (unfractured) floes and add fractured pieces
    for i in range(nfloes2frac, 1, step = -1)
        new_floes = fracture_list[i]
        if !isempty(new_floes)
            n_new_floes = length(new_floes)
            new_floes.id .= range(max_floe_id + 1, max_floe_id + n_new_floes)
            push!.(new_floes.parent_ids, floes.id[frac_idx[i]])
            append!(floes, new_floes)
            max_floe_id += n_new_floes
            StructArrays.foreachfield(f -> deleteat!(f, frac_idx[i]), floes)
        end
    end
    return max_floe_id
end

