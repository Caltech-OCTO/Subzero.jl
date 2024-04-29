

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
    HiblerYieldCurve{FT<:AbstractFloat}<:AbstractFractureCriteria

Type of AbstractFractureCriteria that creates a yield curve that determines if a
floe fractures based off if its stress in principal stress space  is inside or
outside of the yield curve.
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
    HiblerYieldCurve(::Type{FT}, args...)

A float type FT can be provided as the first argument of any HiblerYieldCurve
constructor. A HiblerYieldCurve of type FT will be created by passing all
other arguments to the correct constructor. 
"""
HiblerYieldCurve(::Type{FT}, args...) where {FT <: AbstractFloat}=
    HiblerYieldCurve{FT}(args...)

"""
    HiblerYieldCurve(args...)

If a type isn't specified, HiblerYieldCurve will be of type Float64 and the
correct constructor will be called with all other arguments.
"""
HiblerYieldCurve(args...) = HiblerYieldCurve{Float64}(args...)

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
MohrsCone{FT<:AbstractFloat}<:AbstractFractureCriteria

Type of AbstractFractureCriteria that creates a cone in principal stress space
that determines if a floe fractures based off if its stress in principal stress
space  is inside or outside of the cone.
Fields:
    verticies   <PolyVec> vertices of criteria in principal stress space
Note:
    Concepts from the following papter -
    Weiss, Jérôme, and Erland M. Schulson. "Coulombic faulting from the grain
    scale to the geophysical scale: lessons from ice." Journal of Physics D:
    Applied Physics 42.21 (2009): 214017.
"""
struct MohrsCone{FT<:AbstractFloat}<:AbstractFractureCriteria
    vertices::PolyVec
    
    function MohrsCone{FT}(
        vertices::PolyVec
    ) where {FT <: AbstractFloat}
        try
            valid_polyvec!(vertices)
        catch
            throw(ArgumentError("The given vertices for the Mohr's Cone can't \
            be made into a valid polygon and thus the initial yield \
            curve can't be created."))
        end
        new{FT}(vertices)
    end
end

"""
    MohrsCone(::Type{FT}, args...)

A float type FT can be provided as the first argument of any MohrsCone
constructor. A MohrsCone of type FT will be created by passing all
other arguments to the correct constructor. 
"""
MohrsCone(::Type{FT}, args...) where {FT <: AbstractFloat}=
    MohrsCone{FT}(args...)

"""
    MohrsCone(args...)

If a type isn't specified, MohrsCone will be of type Float64 and the correct
constructor will be called with all other arguments.
"""
MohrsCone(args...) = MohrsCone{Float64}(args...)

"""
    calculate_mohrs(σ1, σ2, σ11, σ22)

Creates PolyVec from vertex values for Mohr's Cone (triangle in 2D)
Inputs:
    σ1  <AbstractFloat> x-coordiante of first point in cone
    σ2  <AbstractFloat> y-coordiante of first point in cone
    σ11 <AbstractFloat> x-coordinate of one vertex of cone and negative of the
            y-coordinate of adjacend vertex in principal stress space
    σ22 <AbstractFloat> y-coordinate of one vertex of cone and negative of the
    x-coordinate of adjacend vertex in principal stress space
Output:
    Mohr's Cone vertices (triangle since we are in 2D) in principal stress space
"""
calculate_mohrs(σ1, σ2, σ11, σ22) = valid_polyvec!([[
    [σ1, σ2],
    [σ11, σ22],
    [σ22, σ11],
    [σ1, σ2],
]])

"""
    calculate_mohrs(
        q,
        σc,
        σ11;
        σ1 = nothing,
        σ2 = nothing,
        σ22 = nothing,
    )

Calculate Mohr's Cone coordinates in principal stress space.
Inputs:
    q   <AbstractFloat> based on the coefficient of internal friction (µi) by
            ((μi^2 + 1)^(1/2) + μi^2
    σc  <AbstractFloat> uniaxial compressive strength
    σ11 <AbstractFloat> negative of the x-coordinate of one vertex of cone
            (triangle in 2D) and negative of the y-coordinate of adjacend vertex
            in principal stress space
Outputs:
    Mohr's Cone vertices (triangle since we are in 2D) in principal stress space
Note:
    Concepts from the following papter -
    Weiss, Jérôme, and Erland M. Schulson. "Coulombic faulting from the grain
    scale to the geophysical scale: lessons from ice." Journal of Physics D:
    Applied Physics 42.21 (2009): 214017.
    Equations taken from original version of Subzero written in MATLAB
"""
function calculate_mohrs(
    q = 5.2,
    σc = 2.5e5,
    σ11 = -3.375e4,
)
    σ1 = ((1/q) + 1) * σc / ((1/q) - q)
    σ2 = q * σ1 + σc
    σ22 = q * σ11 + σc
    return calculate_mohrs(-σ1, -σ2, -σ11, -σ22)
end

"""
    MohrsCone{FT}(val::AbstractFloat, args...)

Calculate Mohr's Cone vertices given calculate_mohrs arguments.
"""
MohrsCone{FT}(args...) where FT = MohrsCone{FT}(calculate_mohrs(args...))

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
    return
end

"""
    update_criteria!(::MohrsCone, floes)

Mohr's cone is not time or floe dependent so it doesn't need to be updates.
"""
function update_criteria!(::MohrsCone, floes)
    return
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
    floe_settings,
    Δt,
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
                floe_settings,
                rng,
            )
            conserve_momentum_change_floe_shape!(
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
        rng,
        fracture_settings,
        floe_settings,
        Δt,
    )
Splits a given floe into pieces using voronoi tesselation.
User will recieve a warning if floe isn't split.
Inputs:
    floe              <Floe> floe in simulation
    rng               <RNG> random number generator used for voronoi tesselation
    fracture_settings <FractureSettings> simulation's fracture settings
    floe_settings     <FloeSettings> simulation's settings for making floes
    Δt                <Int> length of simulation timesteps in seconds
Outputs:
    new_floes   <StructArray{Floes}> list of pieces floe is split into, each of
                    which is a new floe
"""
function split_floe(
    floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    rng,
    fracture_settings,
    floe_settings,
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
        floe_poly = LG.Polygon(rmholes(floe.coords))
        pieces_polys = Vector{LG.Polygon}()
        for p in pieces
            append!(
                pieces_polys,
                get_polygons(
                    LG.intersection(LG.Polygon(p), floe_poly)
                ),
            )
        end
        # Conserve mass within pieces
        pieces_areas = [LG.area(p) for p in pieces_polys]
        total_area = sum(pieces_areas)
        # Create floes out of each piece
        for i in eachindex(pieces_polys)
            if pieces_areas[i] > 0
                mass = floe.mass * (pieces_areas[i]/total_area)
                height = mass / (floe_settings.ρi * pieces_areas[i])
                pieces_floes = poly_to_floes(
                    FT,
                    pieces_polys[i],
                    height,
                    0;  # Δh - range of random height difference between floes
                    floe_settings = floe_settings,
                    rng = rng,
                    u = floe.u,
                    v = floe.v,
                    ξ = floe.ξ,
                )
                append!(new_floes, pieces_floes)
            end
        end
    end

    # Conserve momentum and update strain
    conserve_momentum_fracture_floe!(
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
        floe_settings,
        Δt,
    )
Fractures floes that meet the criteria defined in the fracture settings.
Inputs:
    floes       <StructArray{Floe}> model's list of floes
    max_floe_id <Int> maximum ID of any floe created so far in simulation
    rng         <RNG> random number generator
    fracture_settings   <FractureSettings> sim's fracture settings
    floe_settings       <FloeSettings> sim's settings to make floes
    Δtout               <Int> length of simulation timestep in seconds
Outputs:
    max_floe_id <Int> new highest floe ID after adding new floes to floe array.
    Floe pieces added to floe array and original fractured floes removed.
"""
function fracture_floes!(
    floes::StructArray{<:Floe{FT}},
    max_floe_id,
    rng,
    fracture_settings,
    floe_settings,
    Δt,
) where {FT <: AbstractFloat}
    # Determine which floes will fracture
    frac_idx = determine_fractures(
        floes,
        fracture_settings.criteria,
        floe_settings.min_floe_area,
    )
    # Initialize list for new floes created from fracturing existing floes
    nfloes2frac = length(frac_idx)
    fracture_list = [StructArray{Floe{FT}}(undef, 0) for _ in 1:nfloes2frac]
    # Fracture floes that meet criteria 
    Threads.@threads for i in 1:nfloes2frac
        # Deform floe around largest impact site
        if fracture_settings.deform_on
            fidx = frac_idx[i]
            max_overlap = FT(0)
            max_overlap_idx = 0
            for i in 1:floes.num_inters[fidx]
                if (
                    floes.interactions[fidx][i, floeidx] > 0 &&
                    max_overlap < floes.interactions[fidx][i, overlap]
                )
                    max_overlap = floes.interactions[fidx][i, overlap]
                    max_overlap_idx = i
                end
            end
            if max_overlap > 0
                deforming_inter = floes.interactions[fidx][max_overlap_idx, :]
                deforming_floe_idx = Int(deforming_inter[floeidx])
                if deforming_floe_idx <= length(floes)
                    deform_floe!(
                        LazyRow(floes, frac_idx[i]), 
                        floes.coords[deforming_floe_idx],
                        deforming_inter[xforce:yforce],
                        floe_settings,
                        Δt,
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
            floe_settings,
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

