@kwdef struct Constants{FT<:AbstractFloat}
    ρo::FT = 1027.0             # Ocean density
    ρa::FT = 1.2                # Air density
    Cd_io::FT = 3e-3            # Ice-ocean drag coefficent
    Cd_ia::FT = 1e-3            # Ice-atmosphere drag coefficent
    Cd_ao::FT = 1.25e-3         # Atmosphere-ocean momentum drag coefficient
    f::FT = 1.4e-4              # Ocean coriolis frequency
    turnθ::FT = 15π/180         # Ocean turn angle
    L::FT = 2.93e5              # Latent heat of freezing [Joules/kg]
    k::FT = 2.14                # Thermal conductivity of surface ice[W/(m*K)]
    ν::FT = 0.3                 # Poisson's ratio
    μ::FT = 0.2                 # Coefficent of friction
    E::FT = 6e6                 # Young's Modulus
end

"""
Constants(::Type{FT}, args...)

A float type FT can be provided as the first argument of any Constants
constructor. A Constants of type FT will be created by passing all other
arguments to the correct constructor. 
"""
Constants(::Type{FT} = Float64; kwargs...) where {FT <: AbstractFloat} =
    Constants{FT}(args...; kwargs...)