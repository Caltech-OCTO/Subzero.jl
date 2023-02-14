"""
Structs and functions for fracturing floes
"""

abstract type AbstractFractureCriteria{FT<:AbstractFloat} end

struct HiblerYieldCurve{FT}<:AbstractFractureCriteria{FT}
    pstar::FT
    c::FT
end

struct MohrCone{FT}<:AbstractFractureCriteria{FT}
    sig1::FT
    sig2::FT
    sig3::FT
    sig4::FT
end

struct NoFracture{FT}<:AbstractFractureCriteria{FT} end

