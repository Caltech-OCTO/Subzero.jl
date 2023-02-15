"""
Structs and functions for fracturing floes
"""

abstract type AbstractFractureCriteria{FT<:AbstractFloat} end

struct HiblerYieldCurve{FT}<:AbstractFractureCriteria{FT}
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
    vertices = splitdims([x'; y'])
    vertices = rotate_degrees(vertices, 45)
    vertices = translate(vertices, fill(-p/2, 2))
    return vertices
end

HiblerYieldCurve(floes, pstar = 2.25e5, c = 20.0) =
    HiblerYieldCurve(pstar, c, calculate_hibler(floe, pstar, c))

struct MohrCone{FT}<:AbstractFractureCriteria{FT}
    sig1::FT
    sig2::FT
    sig3::FT
    sig4::FT
    vertices::PolyVec{FT}
end

struct NoFracture{FT}<:AbstractFractureCriteria{FT} end

function update_criteria!(criteria::HiblerYieldCurve, floes)
    criteria.vertices = calculate_hibler(floes, criteria.pstar, criteria.c)
end

function determine_fractures(floes, criteria::AbstractFractureCriteria, min_floe_area = 0)
    # Determine if floe stresses are in or out of criteria allowable regions
    update_criteria!(criteria, floes)
    σ_min, σ_max = extrema.(eigvals.(floes.stress))
    x, y = seperate_xy(criteria.vertices)
    in_on = inpoly2(hcat(σ_min, σ_max), hcat(x, y))
    # If stresses are outside of criteria regions, fracture the floe
    frac_idx = ![in_on[:, 1] .|  in_on[:, 2], :]
    frac_idx[floes.area .< min_floe_area] .= false
    return frac_idx
end

function deform_floe!(floe, deformer_coords, deforming_forces)
    poly = LG.Polygon(floe.coords)
    deformer_poly = LG.Polygon(deformer_coords)
    overlap_region = sortregions(LG.intersection(poly, deformer_poly))[1]
    if LG.area(overlap_region) > 0
        rcent = LG.GeoInterface.coordinates(LG.centroid(overlap_region))::Vector{Float64}
        dist = calc_point_poly_dist(rcent[1], rcent[2], LG.GeoInterface.coordinates(overlap_region))
        force_fracs = deforming_forces ./ 2norm(deforming_forces)
        Δx, Δy = abs(dist) .* force_fracs
        deformer_poly = LG.Polygon(translate(deformer_coords, [Δx, Δy]))
        new_floe = sortregions(LG.difference(poly, deformer_poly))[1]
        new_floe_area = LG.area(new_floe)
        if new_floe_area > 0 && new_floe_area/floe.area > 0.9
            new_floe_centroid = LG.GeoInterface.coordinates(LG.centroid(overlap_region))::Vector{Float64}
            floe.centroid = new_floe_centroid
            floe.coords = LG.GeoInterface.coordinates(new_floe)::PolyVec{Float64}
            floe.area = new_floe_area
        end
    end
end

function split_floe(floe, npieces, rng, ::Type{T} = Float64) where T
    new_floes = StructArray{Floe{T}}(undef, 0)
    floe_poly = LibGEOS.Polygon(floe.coords)
    scale_fac = fill(2floe.rmax, 2)
    trans_vec = [floe.centroid[1] - floe.rmax, floe.centroid[2] - floe.rmax]
    pieces = generate_voronoi_coords(npieces, scale_fac, trans_vec, floe.coords, rng; t = T)
    for p in pieces
        piece_poly = LG.intersection(LG.Polygon(p), floe_poly)
        pieces_floes = poly_to_floes(piece_poly, floe.h, 0; ρi = ρi, u = floe.u, v = floe.v, ξ = floe.ξ,
                                                            mc_n = mc_n, rng = rng, t = T)
        append!(new_floes, pieces_floes)
    end

    new_floes.strain .= floe.strain
    new_floes.p_dxdt .= floe.p_dxdt
    new_floes.p_dydt .= floe.p_dydt
    new_floes.p_dudt .= floe.p_dudt
    new_floes.p_dvdt .= floe.p_dvdt
    new_floes.p_dξdt .= floe.p_dξdt

    return new_floes
end

function fracture_floes(floes, criteria::AbstractFractureCriteria, rng, npieces = 3, min_floe_area = 0, ::Type{T} = Float64) where T
    frac_idx = determine_fractures(floes, min_floe_area, criteria)
    fractured_list = Vector{StructVector{Floes}}(undef, sum(frac_idx))
    for i in frac_idx
        ifloe = floes[i]
        inters = ifloe.interactions
        inters = inters[!isinf(inters[:, floeidx])]
        if !isempty(inters)
            _, max_inters_idx = findmax(inters[:, overlap])
            deforming_inter = inters[max_inters_idx, :]
            deforming_floe_idx = floes[deforming_inter[floeidx]]
            if deforming_floe_idx <= length(floes)
                deform_floe!(ifloe, floes.coords[deforming_floe_idx], deforming_inter[xforce:yforce])
            end
        end
        fractured_list[i] = split_floe(ifloe, npieces, rng, T)
        floes[i] = ifloe
    end
    for idx in frac_idx
        if !isempty(fractured_list[idx])
            StructArrays.foreachfield(f -> deleteat!(f, idx), floes)
            append!(floes, fractured_list[idx])
            # TODO: need to set IDs
        end
    end
end



