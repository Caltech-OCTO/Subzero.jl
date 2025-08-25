export NoFracture, HiblerYieldCurve, MohrsCone, FractureSettings

"""
    abstract type AbstractFractureCriteria

Abstract super type for fracture criteria, which determines if a floe fractures given its
physical condition. Right now, fracture criteria define a set of vertices that define a
shape in principal stress space. The minimum and maximum eigenvalues of a floe's stress field
are then a point within the principal stress space. If that stress point falls outside of
the criteria-verticies defined vertices, it is a stress big enough to fracture the floe.
Otherwise the floe will not be fractured. If another form of fracture criteria is desired,
some additional code may need to be pulled into helper fuctions within the fracture.jl file
and then dispatched off of using subtypes of AbstractFractureCriteria. This is encuarged if It
is helpful!

Right now, the subtypes serve as dispatch types for the following two methods.

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
struct NoFracture<:AbstractFractureCriteria end

#= No floes are fractured when the fracture criteria is `NoFracture` and this function 
won't even be called due to given the `fractures_on` keyword is `false`. It is defined for
API completeness. =#
_determine_fractures(::NoFracture, _, _) = nothing


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
    poly::Polys{FT}
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
    # TODO: eventually make with SVectors! 
    poly = GI.Polygon([ring_coords])
    return _move_poly(FT, poly, -p/2, -p/2,  π/4)
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
    HiblerYieldCurve struct with vertices determined using the _calculate_hibler
    function.
"""
HiblerYieldCurve{FT}(
    floes::StructArray{<:Floe{FT}},
    pstar = 2.25e5,
    c = 20.0,
) where {FT <: AbstractFloat} =
    HiblerYieldCurve{FT}(
        pstar,
        c,
        _calculate_hibler(FT, mean(floes.height), pstar, c),
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
    poly::Polys{FT}
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
    MohrsCone{FT}(val::AbstractFloat, args...)

Calculate Mohr's Cone vertices given _calculate_mohrs arguments.
"""
MohrsCone{FT}(args...) where FT = MohrsCone{FT}(_calculate_mohrs(FT, args...))

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
function update_criteria!(criteria::HiblerYieldCurve{FT}, floes) where FT
    criteria.poly = _calculate_hibler(
        FT,
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
    FractureSettings{CT}

Settings needed for fractures within the model. 

## _Fields_
- `fractures_on::Bool`: if true, floes will fracture, else they will not (Default = false)
- `criteria::CT`: defines which fracture criteria (subtype of [`AbstractFractureCriteria`](@ref)) are used to determine which floes to fracture
- `Δt::Int`: determines how many simulation timesteps between attempting floe fracture (Default = 0)
- `deform_on::Bool`: if true, then the floe will be deformed around floe primarily causing the fracture, identified by the largest overlap area on the most recent set of collisions (Default = False)
- `npieces::Int`: how many pieces to try to split a fractured floe into (Default = 3)

## _Keyword arguments_
- Each of the above fields is an optional keyword argument.

## _Examples_
- Creating default `FractureSettings` 
```jldoctest
julia> FractureSettings()
FractureSettings{NoFracture}
  ⊢ fractures_on = false
  ⊢ criteria = NoFracture
  ⊢ Δt = 0
  ⊢ deform_on = false
  ⊢ npieces = 3
```

- Creating a Float32 `FractureSettings` with fractures on
```jldoctest
julia> mohrs_cone_criteria = MohrsCone(Float32);

julia> FractureSettings(; fractures_on = true, criteria = mohrs_cone_criteria)
FractureSettings{MohrsCone{Float32}}
  ⊢ fractures_on = true
  ⊢ criteria = MohrsCone{Float32}
  ⊢ Δt = 0
  ⊢ deform_on = false
  ⊢ npieces = 3
```
"""
@kwdef struct FractureSettings{CT<:AbstractFractureCriteria}
    fractures_on::Bool = false
    criteria::CT = NoFracture()
    Δt::Int = 0
    deform_on::Bool = false
    npieces::Int = 3

    function FractureSettings{CT}(
        fractures_on,
        criteria::CT,
        Δt,
        deform_on,
        npieces,
    ) where {CT <: AbstractFractureCriteria}
        if fractures_on
            if Δt < 0
                @warn "Fracturing can't occur with negative timesteps. Turning \
                    fracturing off."
                fractures_on = false
            elseif criteria isa NoFracture
                @warn "Fracturing can't occur on with NoFracture criteria. \
                    Turning fracturing off."
                fractures_on = false
            elseif npieces < 2
                @warn "Fracturing can't occur on with npieces < 2 as this \
                    won't split floe. Turning fracturing off."
                fractures_on = false
            end
        end
        if !fractures_on && deform_on
            @warn "Deformation can't occur on without fracturing. Turning \
                deformation off."
            deform_on = false
        end
        new{CT}(fractures_on, criteria, Δt, deform_on, npieces)
    end
    FractureSettings(
        fractures_on,
        criteria::CT,
        Δt,
        deform_on,
        npieces,
    ) where {CT <: AbstractFractureCriteria} = 
        FractureSettings{CT}(
            fractures_on,
            criteria,
            Δt,
            deform_on,
            npieces,
        )
end

# Pretty printing for FractureSettings showing key dimensions
function Base.show(io::IO, fracture_settings::FractureSettings{CT}) where {CT}
    overall_summary = "FractureSettings{$CT}"
    consts_summary = ""
    for name in fieldnames(FractureSettings)
        if name == :criteria
            consts_summary *= "  ⊢" * " " * (string(name) * " = " * string(typeof(getfield(fracture_settings, name))) * "\n")
        else
            consts_summary *= "  ⊢" * " " * (string(name) * " = " * string(getfield(fracture_settings, name)) * "\n")
        end
    end
    print(io, overall_summary, "\n", consts_summary)
end