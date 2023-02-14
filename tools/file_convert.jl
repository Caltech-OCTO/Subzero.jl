"""
Code for converting between .MAT and .jld2 filetypes
"""
# Uses JLD2.jl and SplitApplyCombine.jl

# Given .mat file with a struct that holds an array of polyshapes.
# This doesn't open well with Julia because polyshape doesn't exist. 
# Do a version of the following commented-out code in MATLAB.
# It opens the file and resaves just the vertices in a cell array:

# F = load('FloeShapes.mat');
# polyVertices = cell(size(F, 2), 1)
# for i = 1:size(F.poly, 2)
#    polyVertices(i, 1) = {F.poly(i).Vertices}
# end
# save('FloeVertices.mat', 'polyVertices')

# Then switch to Julia, open the new file, transform the coordiantes into
# PolyVecs using SplitApplyCombine and re-save as a JLD2 file for use
vertices_data = matread("Desktop/FloeVertices.mat")
vertices = [[splitdims(v')] for v in vertices_data["polyVertices"]]
jldsave("Desktop/floe_shapes.jld2"; floe_vertices = vertices)


vars = matread("test/inputs/floe2.mat")
coords = vars["floe"]["c_alpha"]
coords = splitdims(coords)
coords = Subzero.translate([coords], [vars["floe"]["Xi"], vars["floe"]["Yi"]])
floe = Floe(coords, vars["floe"]["h"], 0.0)
floe.centroid = [vars["floe"]["Xi"], vars["floe"]["Yi"]]
floe.mc_x = vec(vars["floe"]["X"])[vec(vars["floe"]["A"]) .== 1]
floe.mc_y = vec(vars["floe"]["Y"])[vec(vars["floe"]["A"]) .== 1]
floe.u =vars["floe"]["Ui"]
floe.v =vars["floe"]["Vi"]
floe.ξ = vars["floe"]["ksi_ice"]
floe.α = vars["floe"]["alpha_i"]
floe.interactions = vars["floe"]["interactions"]
floe.fxOA = vars["floe"]["FxOA"]
floe.fyOA = vars["floe"]["FyOA"]
floe.trqOA = vars["floe"]["torqueOA"]
floe.collision_force = vars["floe"]["collision_force"]
floe.collision_trq = vars["floe"]["collision_torque"]
floe.stress = vars["floe"]["Stress"]
floe.strain = vars["floe"]["strain"]
for x in eachslice(vars["floe"]["StressH"], dims = 3)
    if sum(x) != 0
        push!(floe.stress_history, x)
    end
end

jldsave("test/inputs/test_floes.jld2"; stress_strain_floe = floe)