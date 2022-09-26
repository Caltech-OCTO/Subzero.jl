"""
Main module for `Subzero.jl` - UW's sea ice model ported from MATLAB to Julia
"""
module Subzero

export 
    Model,
    Simulation,
    run!

import LibGEOS as LG
using Plots

# some structs here!

include("floe_operations.jl")
include("model.jl")

end
