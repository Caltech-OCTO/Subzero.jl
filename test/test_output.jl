

# Test checkpoint output writer
file = jldopen("output/sim/checkpoint.jld2", "r")
@test Set(keys(file)) == Set(["ocean", "metadata", "atmos", "floes"])
@test typeof(file["ocean/0"]) == Subzero.Ocean
@test typeof(file["atmos/0"]) == Subzero.Atmos
@test typeof(file["floes/0"]) == StructArray{Subzero.Floes}

