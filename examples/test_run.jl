using JLD2, Random, SplitApplyCombine, Statistics, StructArrays, Subzero
import LibGEOS as LG

frac_stress = [-29955.396 -3428.008; -3428.008	-1942.0464]
frac_deform_floe = Floe(
    [[
        [-50548.0, -49996.0],
        [-50551.00000000001, -37790.0],
        [-20856.0, -32519.0],
        [-20930.0, -49990.0],
        [-50548.0, -49996.0],
    ]],
    0.25,
    0.0,
    u = 0.1,
    v = -0.2,
    ξ = 0.05,
)
frac_floe = deepcopy(frac_deform_floe)  # Without interactions, won't deform
no_frac_floe = Floe(  # This floe is colliding with frac_deform_floe
    [[
        [1467.795, -25319.563],
        [1664.270, -25640.216],
        [-1105.179, -25640.216],
        [-17529.019, -50035.583],
        [-21193.828, -50088.777],
        [-21370.170, -32618.322],
        [-21247.656, -31077.536],
        [-12818.593, -27031.048],
        [1467.795, -25319.563],
    ]],
    0.25,
    0.0,
)
no_frac_small = Floe(  # This floe is too small to fracture or deform
    [[
        [1e3, 1e3],
        [1e3, 1.5e3],
        [1.5e3, 1.5e3],
        [1.5e3, 1e3],
        [1e3, 1e3],
    ]],
    0.25,
    0.0,
)
frac_deform_floe.stress = frac_stress
frac_deform_floe.interactions = collect([
    3,
    -279441968.984,
    -54223517.438,
    -21091.0918258529,
    -40358.0042297616,
    -148920620521.112,
    6795329.38154967,
]')
frac_deform_floe.p_dudt = 0.11
frac_floe.stress = frac_stress
no_frac_small.stress = frac_stress

floes = StructArray([
    frac_deform_floe, frac_floe, no_frac_floe, no_frac_small
])
floes.id .= collect(1:4)
frac_settings = FractureSettings(
    fractures_on = true,
    criteria = HiblerYieldCurve(floes),
    Δt = 75,
    deform_on = true,
)

# Test determine_fractures
frac_idx = Subzero.determine_fractures(
        floes,
        HiblerYieldCurve(floes),
        1e6
)

# Test deform_floe!
floe1_copy = deepcopy(floes[1])
colliding_coords = no_frac_floe.coords
deforming_forces = frac_deform_floe.interactions[xforce:yforce]
init_overlap = LG.area(LG.intersection(
    LG.Polygon(floe1_copy.coords),
    LG.Polygon(colliding_coords),
))
Subzero.deform_floe!(floe1_copy, colliding_coords, deforming_forces)

# Test split_floe
new_floes = Subzero.split_floe(
    floes[1],
    Xoshiro(3),
    frac_settings,
    CouplingSettings(),
    Constants(),
) 
# Test that the pieces all fit within original floe
og_floe_poly = LG.Polygon(floes.coords[1])
new_floes_polys = LG.MultiPolygon(new_floes.coords)

# Test fracture_floes!
max_idx = Subzero.fracture_floes!(
    floes,
    4,  # start with 4 floes
    Xoshiro(3),
    frac_settings,
    CouplingSettings(),
    SimplificationSettings(),
    Constants(),
)

# User Inputs
const type = Float64::DataType
const Lx = 1e5
const Ly = 1e5
const Δgrid = 10000
const hmean = 0.25
const Δh = 0.0
const Δt = 20

# Model instantiation
grid = RegRectilinearGrid(0, Lx, 0, Ly, Δgrid, Δgrid)
ocean = Ocean(grid, 0.0, 0.0, 0.0)
atmos = Atmos(zeros(grid.dims .+ 1), zeros(grid.dims .+ 1), zeros(grid.dims .+ 1))

# Domain creation - boundaries and topography
nboundary = CollisionBoundary(grid, North())
sboundary = CollisionBoundary(grid, South())
eboundary = CollisionBoundary(grid, East())
wboundary = CollisionBoundary(grid, West())

#island = [[[6e4, 4e4], [6e4, 4.5e4], [6.5e4, 4.5e4], [6.5e4, 4e4], [6e4, 4e4]]]
#topo = TopographyElement([[[-9.5e4, 4.5e4], [-9.5e4, 6.5e4], [-6.5e4, 6.5e4],
#                           [-6.5e4, 4.5e4], [-9.5e4, 4.5e4]]])
#topo_arr = StructVector([topo for i in 1:1])

#topo1 = [[[0, 0.0], [0, 1e5], [2e4, 1e5], [3e4, 5e4], [2e4, 0], [0.0, 0.0]]]
#topo2 = [[[8e4, 0], [7e4, 5e4], [8e4, 1e5], [1e5, 1e5], [1e5, 0], [8e4, 0]]]

#topo_arr = StructVector([TopographyElement(t) for t in [island, topo1, topo2]])

domain = Domain(nboundary, sboundary, eboundary, wboundary)

# Floe instantiation
floe_arr = initialize_floe_field([[[[1.25e4, 8e4], [1.25e4, 6e4], [1.75e4, 6e4], [1.75e4, 8e4], [1.25e4, 8e4]]],
    [[[0.7e4, 8e4], [0.7e4, 6e4], [1.3e4, 6e4], [1.3e4, 8e4], [0.7e4, 8e4]]]], domain, 0.5, 0.0, nhistory = 1)
#floe_arr = initialize_floe_field(50, [1], domain, 0.5, 0.0, rng = Xoshiro(1), nhistory = 1)
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup

modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 50,
    verbose = true,
    fracture_settings = FractureSettings(
        fractures_on = true,
        criteria = HiblerYieldCurve(floe_arr),
        Δt = 1,
        npieces = 3,
        nhistory = 1000,
    ),
    coupling_settings = CouplingSettings(
        coupling_on = false
    )
)

# Output setup
initwriter = InitialStateOutputWriter(dir = "output/sim", filename = "initial_state.jld2", overwrite = true)
gridwriter = GridOutputWriter(50, grid, (10, 10), dir = "output/sim", filename = "grid.nc", overwrite = true)
floewriter = FloeOutputWriter(50, dir = "output/sim", filename = "floes.jld2", overwrite = true)
checkpointwriter = CheckpointOutputWriter(50, dir = "output/sim", overwrite = true)

# Run simulation
run!(simulation, [floewriter, initwriter])

Subzero.create_sim_gif("output/sim/floes.jld2", 
                       "output/sim/initial_state.jld2",
                       "output/sim/test.gif")