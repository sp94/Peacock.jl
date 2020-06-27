# How to load a crystal from the Zoo

## Namespaces

The `Peacock.Zoo` module contains functions to generate some of the photonic crystals that I have studied. To avoid polluting the global namespace, these functions are only accessible after calling `using Peacock.Zoo`, or by calling each function as `Peacock.Zoo.name`.


## Unpacking parameters

Each `make_*` function returns a `NamedTuple` containing the geometry, solver, polarisation, and high symmetry ``k``-points of the crystal. For example, the topological photonic crystal first introduced by Wu *et al* 2015 can be generated using `make_wu_topo`.
```julia
# Load photonic topological insulator from Wu et al 2015
@unpack geometry, solver, polarisation, G, K, M = make_wu_topo(fourier_space_cutoff)

# Preview the crystal
plot(geometry)
```

![](../figures/example_zoo_geometry.png)


## Available crystals

```@autodocs
Modules = [Peacock.Zoo]
Private = false
```
