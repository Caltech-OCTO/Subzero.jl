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