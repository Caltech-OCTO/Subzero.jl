```@meta
CurrentModule = Subzero
```

# Full Subzero API documentation

## Grids

```@docs
AbstractRectilinearGrid
RegRectilinearGrid
```

## Directions
```@docs
AbstractDirection
North
South
East
West
```
## Boundaries
```@docs
AbstractBoundary
OpenBoundary
PeriodicBoundary
CollisionBoundary
MovingBoundary
```
## Topography
```@docs
TopographyElement
initialize_topography_field
```
## Domain
```@docs
Domain
```

## Ocean
```@docs
Ocean
```

## Atmosphere
```@docs
Atmos
```

## Floes
```@docs
Floe
AbstractFloeFieldGenerator
CoordinateListFieldGenerator
VoronoiTesselationFieldGenerator
initialize_floe_field
```

## Model
```@docs
Model
```

## Constants
```@docs
Constants
```

## Physical Process Settings
```@docs
FloeSettings
CouplingSettings
CollisionSettings
FractureSettings
AbstractFractureCriteria
NoFracture
HiblerYieldCurve
MohrsCone
SimplificationSettings
RidgeRaftSettings
WeldSettings
```

## Output Writers
```@docs
AbstractOutputWriter
InitialStateOutputWriter
CheckpointOutputWriter
GridOutputWriter
FloeOutputWriter
OutputWriters
SubzeroLogger
```

## Simulations
```@docs
Simulation
run!
restart!
```

## Plotting
```@docs
plot_sim
```

# Developer-Used Methods

## Simulation Methods
```@docs
timestep_sim!
```
## Output Writer Methods
```@docs
write_data!
write_init_state_data!
write_checkpoint_data!
write_floe_data!
write_grid_data!
calc_eulerian_data!
```

## Developer-Used Types
```@docs
CellFloes
CellStresses
TopographyField
```
