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
    pstar       <AbstractFloat>
    c           <AbstractFloat>
    verticies   <PolyVec> vertices of criteria in principal stress space
"""
mutable struct HiblerYieldCurve{FT<:AbstractFloat}<:AbstractFractureCriteria
    pstar::FT
    c::FT
    vertices::PolyVec{FT}
end

function calculate_hibler(floes, pstar, c)
    compactness = 1
    h = mean(floes.height)
    p = pstar*h*exp(-c*(1-compactness))
    t = range(0, 2π, length = 100)
    a = p*sqrt(2)/2
    b = a/2
    x = a*cos.(t)
    y = b*sin.(t)
    vertices = [splitdims([x'; y'])]
    vertices = rotate_degrees(vertices, 45)
    vertices = translate(vertices, fill(-p/2, 2))
    return vertices
end

HiblerYieldCurve(floes, pstar = 2.25e5, c = 20.0) =
    HiblerYieldCurve(pstar, c, calculate_hibler(floes, pstar, c))

function update_criteria!(criteria::HiblerYieldCurve, floes)
    criteria.vertices = calculate_hibler(floes, criteria.pstar, criteria.c)
end

"""
    determine_fractures(
        floes,
        criteria,
        min_floe_area,
    )

Determines which floes will fracture depending on the criteria provided.
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
    x, y = seperate_xy(criteria.vertices)
    in_on = inpoly2(σvals, hcat(x, y))
    # If stresses are outside of criteria regions, we will fracture the floe
    frac_idx = .!(in_on[:, 1] .|  in_on[:, 2])
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
"""
function deform_floe!(
    floe,
    deformer_coords,
    deforming_forces,
)
    poly = LG.Polygon(floe.coords)
    deformer_poly = LG.Polygon(deformer_coords)
    overlap_region = sortregions(LG.intersection(poly, deformer_poly))[1]
    if LG.area(overlap_region) > 0
        rcent = find_poly_centroid(overlap_region)
        dist = calc_point_poly_dist(
            rcent[1],
            rcent[2],
            find_poly_coords(overlap_region),
        )
        force_fracs = deforming_forces ./ 2norm(deforming_forces)
        Δx, Δy = abs(dist) .* force_fracs
        deformer_poly = LG.Polygon(translate(deformer_coords, [Δx, Δy]))
        new_floe = sortregions(LG.difference(poly, deformer_poly))[1]
        new_floe_area = LG.area(new_floe)
        if new_floe_area > 0 && new_floe_area/floe.area > 0.9
            new_floe_centroid = find_poly_centroid(new_floe)
            floe.centroid = new_floe_centroid
            floe.coords = find_poly_coords(new_floe)
            floe.area = new_floe_area
        end
    end
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
    floe,
    rng,
    fracture_settings,
    coupling_settings,
    consts,
    ::Type{T} = Float64,
) where T
    new_floes = StructArray{Floe{T}}(undef, 0)
    # Generate voronoi tesselation in floe's bounding box
    scale_fac = fill(2floe.rmax, 2)
    trans_vec = [floe.centroid[1] - floe.rmax, floe.centroid[2] - floe.rmax]
    pieces = generate_voronoi_coords(
        fracture_settings.npieces,
        scale_fac,
        trans_vec,
        floe.coords,
        rng;
        t = T
    )
    # Intersect voronoi tesselation pieces with floe
    floe_poly = LG.Polygon(floe.coords)
    for p in pieces
        piece_poly = LG.intersection(LG.Polygon(p), floe_poly)
        pieces_floes = poly_to_floes(
            piece_poly,
            floe.height,
            0;  # Δh - range of random height difference between floes
            ρi = consts.ρi,
            u = floe.u,
            v = floe.v,
            ξ = floe.ξ,
            mc_n = coupling_settings.mc_n,
            nhistory = fracture_settings.nhistory,
            rng = rng,
            t = T,
        )
        append!(new_floes, pieces_floes)
    end

    # Update new floe pieces with parent information
    new_floes.p_dxdt .= floe.p_dxdt
    new_floes.p_dydt .= floe.p_dydt
    new_floes.p_dudt .= floe.p_dudt
    new_floes.p_dvdt .= floe.p_dvdt
    new_floes.p_dξdt .= floe.p_dξdt
    for s in new_floes.strain
        s .= floe.strain
    end

    return new_floes
end

"""
    fracture_floes!(
        floes,
        rng,
        fracture_settings,
        simp_settings,
        max_floe_id,
        ::Type{T} = Float64
    ) where T
Fractures floes that meet the criteria defined in the fracture settings.
Inputs:
    floes   <StructArray{Floe}> model's list of floes
    rng     <RNG> random number generator
    fracture_settings   <FractureSettings> sim's fracture settings
    simp_settings       <SimplificationSettings> sim's simplification settings
    max_floe_id         <Int> highest ID of any floe created in the simulation
Outputs:
    max_floe_id <Int> new highest floe ID after adding new floes to floe array.
    Floe pieces added to floe array and original fractured floes removed.
"""
function fracture_floes!(
    floes,
    max_floe_id,
    rng,
    fracture_settings,
    coupling_settings,
    simp_settings,
    consts,
    ::Type{T} = Float64
) where T
    # Determine which floes will fracture
    frac_idx = determine_fractures(
        floes,
        fracture_settings.criteria,
        simp_settings.min_floe_area,
    )
    nfloes2frac = length(frac_idx)
    # Initialize list for new floes created from fracturing existing floes
    fractured_list = Vector{StructArray{Floe{T}}}(undef, 0)
    # Fracture floes that meet criteria 
    for i in frac_idx
        ifloe = floes[i]
        # Deform floe around largest impact site
        if fracture_settings.deform_on
            inters = ifloe.interactions
            inters = inters[!isinf(inters[:, floeidx])]
            if !isempty(inters)
                _, max_inters_idx = findmax(inters[:, overlap])
                deforming_inter = inters[max_inters_idx, :]
                deforming_floe_idx = deforming_inter[floeidx]
                if deforming_floe_idx <= length(floes)
                    deform_floe!(
                        ifloe, 
                        floes.coords[deforming_floe_idx],
                        deforming_inter[xforce:yforce],
                    )
                end
            end
        end
        # Split flie into pieces
        new_floes = split_floe(
            ifloe,
            rng,
            fracture_settings,
            coupling_settings,
            consts,
            T,
        )
        push!(fractured_list, new_floes)
    end
    # Remove old (unfractured) floes and add fractured pieces
    for i in range(1, nfloes2frac)
        new_floes = fractured_list[i]
        if !isempty(new_floes)
            n_new_floes = length(new_floes)
            new_floes.id .= range(max_floe_id + 1, max_floe_id + n_new_floes)
            new_floes.fracture_id .= floes.id[frac_idx[i]]
            append!(floes, new_floes)
            max_floe_id += n_new_floes
        end
    end
    while !isempty(frac_idx)
        idx = pop!(frac_idx)
        StructArrays.foreachfield(f -> deleteat!(f, idx), floes)
    end
    return max_floe_id
end

