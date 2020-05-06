## Introduction

`Peacock.jl` is a Julia package for studying photonic crystals.
* To do - explain why Peacock
* To do - list features here


## Example usage

Here we reproduce figures 2 and 3 from chapter 5 of _Photonic crystals: molding the flow of light_ by Joannopoulos (see further reading).

### Installation

The package can be installed using the built-in Julia package manager
```julia
using Pkg
Pkg.add("Peacock")
```

After installation, the package can be loaded in Julia
```julia
using Peacock
```


### Defining a photonic crystal in real space

A photonic crystal is defined using the `Geometry` constructor, which requires the optical properties of the crystal, lattice vectors, and resolution
```julia
# Permittivity
function epf(x,y)
    if x^2+y^2 < 0.2^2
        return 8.9
    else
        return 1
    end
end

# Permeability
function muf(x,y)
    return 1
end
    
# Lattice vectors
a1 = [1, 0]
a2 = [0, 1]

# Real-space resolution
d1 = d2 = 0.01

geometry = Geometry(epf, muf, a1, a2, d1, d2)
plot(geometry)
```
> ![Plot of the example geometry](figures/example_plot_geometry.png)

### Defining the solver

The `Solver` approximates the geometry using a truncated Plane Wave Expansion (see further reading). The number of plane waves is determined by the cutoff. Increasing the cutoff will increase the accuracy of the solution, but low-contrast photonic crystals can be well approximated with a relatively low number of plane waves.

```julia
fourier_space_cutoff = 7
solver = Solver(geometry, fourier_space_cutoff)
plot(solver)
```
> ![Plot of the example solver](figures/example_plot_solver_cutoff=7.png)


### Plotting the band structure of a photonic crystal

* Plotting band structure
    * Defining Brillouin zone coordinates

* Visualising the modes


## Further reading

Photonic crystals - Joannaopoulos
http://ab-initio.mit.edu/book/

Methodology - Rumpf CEM
