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

### Model
```@docs
Model
```

## Developer-Used Types
```@docs
CellFloes
CellStresses
TopographyField
```
