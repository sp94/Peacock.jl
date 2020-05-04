## Introduction

`PlaneWaveExpansionMethod.jl` is a Julia package for studying photonic crystals.
* To do - list features here


## Example usage

Here we demonstrate how to reproduce figures 2 and 3 from chapter 5 of _Photonic crystals: molding the flow of light_ by Joannopoulos (see further reading).

### Installation

The package can be installed using the built-in Julia package manager
```julia
using Pkg
Pkg.add("PlaneWaveExpansionMethod")
```

After installation, the package can be loaded in Julia
```julia
using PlaneWaveExpansionMethod
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
```

We can visualise the crystal using the `plot_field` function
```julia
plot_field(g.ep, g.a1, g.a2)
```


### Defining a photonic crystal in Fourier space

The Plane Wave Expansion Method efficiently solves Maxwell's equations for periodic media using a basis of plane waves in Fourier space (see further reading).

We can construct a plane wave basis using the constructor
```julia
plane_wave_cutoff = 7
plane_wave_basis = PlaneWaveBasis(plane_wave_cutoff)
```
The larger the plane wave cutoff the more accurate the simulation will be.

The `Geometry` constructor takes the real-space representation and returns the 
```julia
fourier_space_geometry = FourierSpaceGeometry(geometry, plane_wave_basis)
```


### Plotting the band structure of a photonic crystal

* Plotting band structure
    * Defining Brillouin zone coordinates

* Visualising the modes


## Further reading

Photonic crystals - Joannaopoulos
http://ab-initio.mit.edu/book/

Methodology - Rumpf CEM
