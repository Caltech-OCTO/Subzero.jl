# # HiblerCurveFractureCriteria and MohrsConeFractureCriteria

export AbstractFractureCriteria, NoFracture, HiblerCurveFractureCriteria, MohrsConeFractureCriteria

#= 
## What is fracture criteria? 

Fracture criteria determines if a floe fractures given it's current physical condition.

The `AbstractFractureCriteria` is the abstract type, with `HiblerCurveFractureCriteria`,
`MohrsConeFractureCriteria`, and `NoFracture` as implemented concrete types.

Each of these concrete types represent different method of defining the fracture criteria
(including a lack of criteria for when fractures aren't enabled), and updating the criteria
throughout the run if needed given changing ice pack conditions. A `FractureSettings` struct
requires a concrete subtype of `AbstractFractureCriteria` for the `criteria` field. The
following functions need to be extended for any concrete subtype of
`AbstractFractureCriteria`: `_update_critera!` and `_determine_fractures`.
=#

const FRACTURE_CRITERIA_DEF = ": criteria to determine if a floe fractures given its physical properties."

"""
    abstract type AbstractFractureCriteria

Abstract super type for fracture criteria, which determines if a floe fractures given its
physical condition. The subtypes serve as dispatch types for the following two methods.

## API
The following methods must be implemented for all subtypes:
- `_update_critera!(criteria::AbstractFractureCriteria, floes::StructArray{<:Floe})`
- `_determine_fractures(...)`

`_update_critera!` is called in the ...

`_determine_fractures` is called in the ...

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

# TODO: Add a determine_fractures function that just defaults to nothing

"""
    HiblerCurveFractureCriteria{FT<:AbstractFloat}<:AbstractFractureCriteria

Type of AbstractFractureCriteria that represents a elliptical yield curve in prinicpal
stress space. When checking for floe fracture, the floe's accumulated stress is translated
into a point in principal stress space and then if that point is outside of the
ellipse defined by the `HiblerCurveFractureCriteria` the floe breaks.

A floe's accumulated stress field, `accum_stress`, is translated into principal stress space
by taking the minimum and maximum eigenvalues as its coordinates in principal stress space.

## Fields:
- `pstar::AbstractFloat`: used to tune ellipse for optimal fracturing
- `c::AbstractFloat`: used to tune ellipse for optimal fracturing
- `a::AbstractFloat`: 
- `b::AbstractFloat`:

## Note:
Hibler's 1979 paper as "A Dynamic Thermodynamic Sea Ice Model" says that:
Both `pstar` and `c` relate the ice strength to the ice thickness and ice compactness.
`c` is determined in that 10% open water reduces the strength substantially and `pstar` is
considered a free parameter
"""
@kwdef mutable struct HiblerCurveFractureCriteria{FT<:AbstractFloat}<:AbstractFractureCriteria
    pstar::FT = 2.25e5
    c::FT = 20.0
    poly::Polys{FT}  # TODO: make this defined by a and b
end

"""
    HiblerCurveFractureCriteria(::Type{FT}, args...)

A float type FT can be provided as the first argument of any HiblerCurveFractureCriteria
constructor. A HiblerCurveFractureCriteria of type FT will be created by passing all
other arguments to the correct constructor. 
"""
HiblerCurveFractureCriteria(::Type{FT}, args...; kwargs...) where {FT <: AbstractFloat} =
    HiblerCurveFractureCriteria{FT}(args...; kwargs...)

"""
    HiblerCurveFractureCriteria(args...)

If a type isn't specified, HiblerCurveFractureCriteria will be of type Float64 and the
correct constructor will be called with all other arguments.
"""
HiblerCurveFractureCriteria(args...; kwargs...) = HiblerCurveFractureCriteria{Float64}(args...; kwargs...)

"""
    _calculate_hibler(FT, floes, pstar, c)

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
function _calculate_hibler(::Type{FT}, mean_height, pstar, c) where FT
    compactness = 1  # Could be a user input with future development
    p = pstar*mean_height*exp(-c*(1-compactness))
    α_range = range(zero(FT), FT(2π), length = 100)
    a = p*sqrt(2)/2
    b = a/2
    ring_coords = [(a*cos(α), b*sin(α)) for α in α_range]
    ring_coords[end] = ring_coords[1] # make sure first and last element are exactly the same
    # TODO: just have everything be defined by `a` and `b`!!
    poly = GI.Polygon([ring_coords])
    return _move_poly(FT, poly, -p/2, -p/2,  π/4)
end

"""
    HiblerCurveFractureCriteria(floes, pstar = 2.25e5, c = 20.0)

Calculates Hibler's Elliptical Yield curve using parameters pstar, c, and the
current floe field. 
Inputs:
    floes   <StructArray{Floes}> model's list of floes
    pstar   <AbstractFloat> used to tune ellipse for optimal fracturing
    c       <AbstractFloat> used to tune ellipse for optimal fracturing
Outputs:
    HiblerCurveFractureCriteria struct with vertices determined using the _calculate_hibler
    function.
"""
HiblerCurveFractureCriteria{FT}(
    floes;
    pstar = 2.25e5,
    c = 20.0,
) where {FT <: AbstractFloat} =
    HiblerCurveFractureCriteria{FT}(;
        pstar,
        c,
        poly = _calculate_hibler(FT, mean(floes.height), pstar, c),
    )

"""
    _update_criteria!(criteria::HiblerCurveFractureCriteria, floes)

Update the Hibler Yield Curve vertices based on the current set of floes. The
criteria changes based off of the average height of the floes.
Inputs:
    criteria    <HiblerCurveFractureCriteria> simulation's fracture criteria
    floes       <StructArray{Floe}> model's list of floes
Outputs:
    None. Updates the criteria's vertices field to update new criteria. 
"""
function _update_criteria!(criteria::HiblerCurveFractureCriteria{FT}, floes) where FT
    criteria.poly = _calculate_hibler(
        FT,
        mean(floes.height),
        criteria.pstar,
        criteria.c
    )
    return
end

"""
MohrsConeFractureCriteria{FT<:AbstractFloat}<:AbstractFractureCriteria

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
struct MohrsConeFractureCriteria{FT<:AbstractFloat}<:AbstractFractureCriteria
    poly::Polys{FT}
end

"""
    MohrsConeFractureCriteria(::Type{FT}, args...)

A float type FT can be provided as the first argument of any MohrsConeFractureCriteria
constructor. A MohrsConeFractureCriteria of type FT will be created by passing all
other arguments to the correct constructor. 
"""
MohrsConeFractureCriteria(::Type{FT}, args...) where {FT <: AbstractFloat}=
    MohrsConeFractureCriteria{FT}(args...)

"""
    MohrsConeFractureCriteria(args...)

If a type isn't specified, MohrsConeFractureCriteria will be of type Float64 and the correct
constructor will be called with all other arguments.
"""
MohrsConeFractureCriteria(args...) = MohrsConeFractureCriteria{Float64}(args...)

"""
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
"""
function _calculate_mohrs(::Type{FT}, σ1, σ2, σ11, σ22) where FT
    # TODO: eventually make with SVectors! 
    points = [(σ1, σ2),  (σ11, σ22), (σ22, σ11), (σ1, σ2)]
    return GI.Polygon([points])
end

"""
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
"""
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