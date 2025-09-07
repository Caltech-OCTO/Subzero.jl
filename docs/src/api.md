```@meta
CurrentModule = Subzero
```

# Full Subzero API documentation

!!! warning
    This page is still very much WIP! The documentation, and to some extent the source code, is being cleaned up. This means that right now, some of the documentation is here, and some is in the [tutorial](https://caltech-octo.github.io/Subzero.jl/dev/tutorial/) and the [documentation.md](https://github.com/Caltech-OCTO/Subzero.jl/blob/main/documentation.md) section of the GitHub sections of the documentation website.

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
CheckpointOutputWriter
GridOutputWriter
FloeOutputWriter
InitialStateOutputWriter
OutputWriters
SubzeroLogger
```

## Simulations
```@docs
Simulation
run!
restart!
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
