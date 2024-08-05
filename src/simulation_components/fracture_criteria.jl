# # HiblerCurveFractureCriteria and MohrsConeFractureCriteria

export AbstractFractureCriteria, NoFracture, HiblerCurveFractureCriteria, MohrsConeFractureCriteria

#= 
## What is fracture criteria? 

Fracture criteria determines if a floe fractures given it's current physical condition.

The `AbstractFractureCriteria` is the abstract type, with `HiblerCurveFractureCriteria`,
`MohrsConeFractureCriteria`, and `NoFracture` as currently implemented concrete types.

Each of these concrete types represent different method of defining the fracture criteria
(including a lack of criteria for when fractures aren't enabled), and updating the criteria
throughout the run if needed given changing ice pack conditions. A `FractureSettings` struct
requires a concrete subtype of `AbstractFractureCriteria` for the `criteria` field. The
following functions need to be extended for any concrete subtype of
`AbstractFractureCriteria`: `_update_criteria!` and `_determine_fractures`.
=#

const FRACTURE_CRITERIA_DEF = ": criteria to determine if a floe fractures given its physical properties."

"""
    abstract type AbstractFractureCriteria

Abstract super type for fracture criteria, which determines if a floe fractures given its
physical condition. The subtypes serve as dispatch types for the following two methods.

## API
The following methods must be implemented for all subtypes:
- `_update_criteria!(criteria::AbstractFractureCriteria, floes::StructArray{<:Floe})`
- `_determine_fractures(criteria::AbstractFractureCriteria, floes::StructArray{<:Floe}, floe_settings::FloeSettings)`

`_update_criteria!` is called in the `fracture_floes!` function and takes the fracture
criteria and the current ice floe pack and updates the fracture criteria given the state
of the ice floe pack if needed.

`_determine_fractures` is called in the `fracture_floes!` and takes in the fracture
criteria, the current ice floe pack, and the floe settings and determines which of the ice
floes fracture. It returns a `Boolean` vector equal in length to the `floes` list. If the
`ith` element of the returned vector is `true` then the `i`th floe in `floes` will fracture.
"""
abstract type AbstractFractureCriteria end


#= Default function to NOT update the fracture criteria before determining floe fractures as
the simulation progresses. Any criteria that should not update as the simulation runs can
simply skip implementing this function and fall back on this default. =#
_update_criteria!(::AbstractFractureCriteria, _) = nothing

"""
    NoFracture<:AbstractFractureCriteria

Type of AbstractFractureCriteria representing when fracturing functionality is turned off.
If this is the type provided to the simulation's `FractureSettings`, then fractures will not
occur. As fracturing will not occur, this subtype depends on the default API functions.
"""
struct NoFracture <: AbstractFractureCriteria end

#= No floes are fractured when the fracture criteria is `NoFracture` and this function 
won't even be called due to given the `fractures_on` keyword is `false`. It is defined for
API completeness. =#
_determine_fractures(::NoFracture, _, _) = nothing

abstract type AbstractPrincipalStressScaler end

struct NoPrincipalStressScaler end

_scale_principal_stress!(::NoPrincipalStressScaler, σvals, floe, floe_settings) = return

struct AreaFracStressScaler{FT}
    α::FT
end

function _scale_principal_stress!(stress_scaler::AreaFracStressScaler{FT}, σvals, floe, floe_settings) where FT
    M = (floe.area/floe_settings.min_floe_area).^stress_scaler.α
    σvals .*= M
    return 
end

# Find floe's accumulated stress in principal stress space so that is can be compared to the
# fracture criteria. 
function find_σpoint(floe, floe_settings)
    σvals = eigvals(floe.stress_accum)
    _scale_principal_stress!(floe_settings.stress_calculator, σvals, floe, floe_settings)
    return σvals
end


# - `_scale_principal_stress!(stress_calculator::AbstractStressCalculator{FT}, σvals::Matrix{FT}, floe::FloeType{FT}, floe_settings::FloeSettings)`

#= 
`_scale_principal_stress!` is called within the `find_σpoint` function which is called
within the `_determine_fractures` function. This function takes the stress calculator, the
`floe`'s accumulated stress in prinicpal stress space (`σvals = eigvals(stress_accum)`), the
`floe` itself, and the `floe_settings` and scales `σvals` by some values/ratio using
physical properties within `floe` and `floe_settings`. This is done to approximate changing
the polygon defining the floe's fracture criteria without having to redefine the fracture 
criteria for each floe. This is almost like a proxy for damage. 

    #= Default function to NOT scale stress before checking if stress point in principal stress
space is inside or outside of the fracture criteria polygon. Any calculators that do NOT
want to scale stress in this way can simply skip implementing this function and fall back
on this default. =#
_scale_principal_stress!(::AbstractStressCalculator, σvals, floe, floe_settings) = return

The area-scaling part of the stress calculations comes into play when deciding if a floe
will be fractured. The accumulated stress can be scaled by a ratio of
(floe.area/min_floe_area)^α. By changing the α value, either larger or smaller floes
fracture more eassily. The scaled value will not be saved as this is equivalent to scaling
the simulation fracture criteria (morh's cone or hibler's ellipse, etc.), but it is less
computationally intensive. By leaving the default `α = 0`, this extra scaling will not
take place.
=#

"""
    HiblerCurveFractureCriteria{FT<:AbstractFloat}<:AbstractFractureCriteria

Type of AbstractFractureCriteria that represents a elliptical yield curve in prinicpal
stress space as defined by Hibler in his 1979 paper as "A Dynamic Thermodynamic Sea Ice
Model". When checking for floe fracture, the floe's accumulated stress is translated into a
point in principal stress space and then if that point is outside of the ellipse defined by
the `HiblerCurveFractureCriteria` the floe breaks.

A floe's accumulated stress field, `accum_stress`, is translated into principal stress space
by taking the minimum and maximum eigenvalues as its coordinates in principal stress space.

## Fields:
- `pstar::AbstractFloat`: used to tune ellipse for optimal fracturing (see Note)
- `c::AbstractFloat`: used to tune ellipse for optimal fracturing (see Note)
- `a::AbstractFloat`: 
- `b::AbstractFloat`:
- `h::AbstractFloat`:

The equation of all points within and on the edge of an arbitrary ellipse with the major and
minor axis defined by `a` and `b`, centered on `(h, h)`, and rotated by π/4 radians) is:

ADD LATEXED EQUATION
`((x-h)cos(A) + (y-k)sin(A))^2 / a^2 + ((x-h)sin(A) + (y-k)cos(A))^2/b^2 ≤ 1`

`pstar` and `c` are used to define calculate the values of `a`, `b`, and `h`. Within
Hibler's paper `h = k` and `A = π/4`. Furthermore `a = `

## Note:
Hibler's 1979 paper as "A Dynamic Thermodynamic Sea Ice Model" says that:
Both `pstar` and `c` relate the ice strength to the ice thickness and ice compactness.
`c` is determined in that 10% open water reduces the strength substantially and `pstar` is
considered a free parameter
"""
mutable struct HiblerCurveFractureCriteria{FT<:AbstractFloat, ST<:AbstractPrincipalStressScaler}<:AbstractFractureCriteria
    pstar::FT
    c::FT
    compactness::FT
    a::FT
    b::FT
    h::FT
    stress_scaler::ST
end

const HIBLER_ARG_WARNING = "Creation of a HiblerCurveFractureCriteria requires the user to provide a floe field, a mean ice pack height, or values for a and b to define the elliptical yield curve."

"""
    HiblerCurveFractureCriteria([::Type{FT} = Float64];
        pstar = 2.25e5,
        c = 20.0,
        compactness = 1.0,
        hmean = nothing,
        a = zero(FT),
        b = zero(FT),
        h = zero(FT),
    )

Create a `HiblerCurveFractureCriteria` ... [Add writing...]

"""
function HiblerCurveFractureCriteria(::Type{FT} = Float64;
    pstar = 2.25e5, c = 20.0, compactness = 1.0, floes = nothing,
    hmean = nothing, a = zero(FT), b = zero(FT), h = zero(FT),
    stress_scaler = NoPrincipalStressScaler(),
) where FT
    # If floes input, find mean height of floes
    if !isnothing(floes)
        hmean = mean(floes.height)
    end
    # If either floes or a mean height are provided, calculate `a`, `b`, and `h`
    if !isnothing(hmean)
        a, b, h = _calculate_hibler(FT, hmean, pstar, c, compactness)
    elseif a ==0 || b == 0  # `a` and `b` must be provided if floes or mean height are not
        @warn HIBLER_ARG_WARNING
    end
    return HiblerCurveFractureCriteria{FT}(pstar, c, compactness, a, b, h, stress_scaler)
end

#= Calculate Hibler's Elliptical Yield Curve as described in his 1979 paper "A Dynamic
Thermodynamic Sea Ice Model" given an ice floe pack and a value for `pstar` and `c` as
described above in the HiblerCurveFractureCriteria documentation. =#
function _calculate_hibler(::Type{FT}, hmean, pstar, c, compactness) where FT
    p = pstar * hmean * exp(-c * (1 - compactness))
    a = p * sqrt(2)/2
    b = a/2
    h = -p/2
    return a, b, h
end

#= Update the Hibler Yield Curve vertices based on the current ice floe pack. The criteria
changes based off of the average height of the floes. =#
function _update_criteria!(criteria::HiblerCurveFractureCriteria{FT}, floes) where FT
    criteria.a, criteria.b, criteria.h = _calculate_hibler(FT, mean(floes.height),
        criteria.pstar, criteria.c, criteria.compactness)
    return
end

function _determine_fractures(criteria::HiblerCurveFractureCriteria, floes, floe_settings)
    nfloes = length(floes)
    frac_idx = falses(nfloes)  # TODO: make with an iterator to reduce allocations!
    for i in 1:nfloes
        floe = get_floe(floes, i)
        # if floe is too small, don't fracture
        floe.area ≤ floe_settings.min_floe_area && continue
        # check if principal stress makes it eligible for fracture
        σpoint = find_σpoint(floe, floe_settings)
        if _out_ellipse(σpoint, criteria)
            frac_idx[i] = true
        end
    end
    return range(1, length(floes))[frac_idx]
end

function _out_ellipse((x, y), criteria)
    Δx, Δy = x - criteria.h, y - criteria.h
    pt_val = 0.5 * (((Δx + Δy) / criteria.a)^2 + ((Δx - Δy) / criteria.b)^2)
    return pt_val > 1
end

"""
    MohrsConeFractureCriteria{FT<:AbstractFloat}<:AbstractFractureCriteria

Type of AbstractFractureCriteria that represents a triangle (2D cone) in principal stress
space. When checking for floe fracture, the floe's accumulated stress is translated into a
point in principal stress space and then if that point is outside of the cone defined by
the `MohrsConeFractureCriteria` the floe fractures.

## Fields:
- `poly::Polys{FT}`: polygon of triangle criteria in principal stress space
## Note:
Concepts from the following papter:
Weiss, Jérôme, and Erland M. Schulson. "Coulombic faulting from the grain
scale to the geophysical scale: lessons from ice." Journal of Physics D:
Applied Physics 42.21 (2009): 214017.
"""
struct MohrsConeFractureCriteria{FT<:AbstractFloat, ST<:AbstractPrincipalStressScaler}<:AbstractFractureCriteria
    poly::Polys{FT}
    stress_scaler::ST
end

"""
    MohrsConeFractureCriteria(::Type{FT}, args...)

A float type FT can be provided as the first argument of any MohrsConeFractureCriteria
constructor. A MohrsConeFractureCriteria of type FT will be created by passing all
other arguments to the correct constructor. 
"""
MohrsConeFractureCriteria(::Type{FT}, args...) where {FT <: AbstractFloat} =
    MohrsConeFractureCriteria{FT}(args...)

"""
    MohrsConeFractureCriteria(args...)

If a type isn't specified, MohrsConeFractureCriteria will be of type Float64 and the correct
constructor will be called with all other arguments.
"""
MohrsConeFractureCriteria(args...) = MohrsConeFractureCriteria{Float64}(args...)

#=
    _calculate_mohrs(FT, σ1, σ2, σ11, σ22)

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
=#
function _calculate_mohrs(::Type{FT}, σ1, σ2, σ11, σ22) where FT
    # TODO: eventually make with SVectors! 
    points = [(σ1, σ2),  (σ11, σ22), (σ22, σ11), (σ1, σ2)]
    return GI.Polygon([points])
end

#=
    _calculate_mohrs(
        FT,
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
=#
function _calculate_mohrs(
    ::Type{FT},
    q = 5.2,
    σc = 2.5e5,
    σ11 = -3.375e4,
) where FT
    σ1 = ((1/q) + 1) * σc / ((1/q) - q)
    σ2 = q * σ1 + σc
    σ22 = q * σ11 + σc
    return _calculate_mohrs(FT, -σ1, -σ2, -σ11, -σ22)
end

"""
    MohrsConeFractureCriteria{FT}(val::AbstractFloat, args...)

Calculate Mohr's Cone vertices given _calculate_mohrs arguments.
"""
MohrsConeFractureCriteria{FT}(args...) where FT = MohrsConeFractureCriteria{FT}(_calculate_mohrs(FT, args...))


"""
    _determine_fractures(
        floes,
        criteria,
        min_floe_area,
    )

Determines which floes will fracture depending on the principal stress criteria.
Inputs:
    floes           <StructArray{Floe}> model's list of floes
    criteria        <AbstractFractureCriteria> fracture criteria
    floe_settings   <FloeSettings> Floe settings. Contains Floe properties and stress 
                    calculator.
Outputs:
    <Vector{Int}> list of indices of floes to fracture 
"""
function _determine_fractures(criteria::MohrsConeFractureCriteria, floes, floe_settings)
    # If stresses are outside of criteria regions, we will fracture the floe
    frac_idx = [!GO.coveredby(find_σpoint(get_floe(floes, i), floe_settings), criteria.poly) for i in eachindex(floes)]
    frac_idx[floes.area .< floe_settings.min_floe_area] .= false
    return range(1, length(floes))[frac_idx]
end


#=
TODO: 
- Switich Mohr's cone to use triangle endpoints
- rewrite _determine_fractures to take advantage of above simplifications and have it dispatch
- Create an "area-scaled" Hiber's ellipse fracture criteria
- Remove the principal stress scaling from the stress calculators
- Clean up file
- Move tests relating to the fracture critera to the test fracture criteria file
=#