<div align="center">
<img src="logo.png" alt="Peacock logo"></img>
</div>

## Introduction

`Peacock.jl` - or the **P**lane-wave **E**xpansion **A**pproach to **C**haracterising **O**ptical **C**rystals in **k**-space - is a Julia package for studying photonic crystals using the Plane Wave Expansion Method [references].

Photonic crystals are materials whose optical properties arise from the structuring of the material when the size of the structures are comparable to the wavelengths of light. `Peacock.jl` is named for the irridescent colours of peacock feathers which arise not from pigmentation but from their photonic crystal structure, as shown below.

![Different zooms of a Peacock](zoom.png)
<sub>Credits: Peacock by allanlau2000 from pixabay. Feather by suju from pixabay. Electron microscope image of photonic crystal structure from Zi, Jian, et al. "Coloration strategies in peacock feathers.",  *Proceedings of the National Academy of Sciences* 100.22 (2003): 12576-12578. Copyright (2003) National Academy of Sciences.</sub>

As well as occuring naturally as in animals such as peacocks, advances in nanofabrication mean that 'designer' photonic crystals can now be manufactured for unprecedented control over the flow of light, with applications ranging from optical fibers [refs] to photonic circuitry [refs]. Photonic crystals are also a promising platform for more exotic materials like topological insulators [refs].

* Solve for...
	* Transverse electric (TE) and transverse magnetic (TM) modes of 2D photonic crystals
	* Non-orthogonal unit cells
	* Inhomogeneous permittivity and/or permeability
	* Chern numbers of topological photonic crystals using built-in Wilson loop methods

* Focused on ease of use
	* Install with one line in Julia's package manager
	* Simple visualisation of geometry, fields, and fully labelled band diagrams



## Example usage

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

A photonic crystal is defined using the `Geometry` constructor, which requires the optical properties of the crystal, lattice vectors, and resolution. Here we demonstrate how to define the photonic crystal used in figures 2 and 3 of chapter 5 of _Photonic crystals: molding the flow of light_ by Joannopoulos (see further reading).
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

Band structures can be generated using `plot_band_diagram(solver, k_path, polarisation)`. If the elements of `k_path` are instances of `BrillouinZoneCoordinate` then the band diagram will be fully labelled.

```julia
G = BrillouinZoneCoordinate(  0,   0, "Γ")
X = BrillouinZoneCoordinate(1/2,   0, "X")
M = BrillouinZoneCoordinate(1/2, 1/2, "M")
k_path = [G,X,M,G]

figure(figsize=(4,3))
plot_band_diagram(solver, k_path, TE, color="red",
            bands=1:4, dk=0.1, frequency_scale=1/2pi)
plot_band_diagram(solver, k_path, TM, color="blue",
            bands=1:4, dk=0.1, frequency_scale=1/2pi)
ylim(0,0.8)
```
> ![Plot of a band structure](figures/example_plot_band_diagram.png)

This reproduces figure 2 of chapter 5 of Joannopoulos *et al* 2008.

### Plotting the modes

`Mode`s can be visualised using `plot(mode)`. By default the full Bloch wave is plotted - set `bloch_phase=false` to plot the cell-periodic part of Bloch mode.

The out of plane field component is plotted - for TE and TM polarisations this will be the magnetic and electric fields, respectively. The titles of the figures are set automatically using the `label` of the `Mode` object.

```julia
frequencies, modes = solve(solver, X, TM)
plot(modes[2], bloch_phase=true)
plot(modes[2], bloch_phase=false)
```

> ![Plot of some modes](figures/example_plot_modes.png)

This reproduces figure 3 of chapter 5 of Joannopoulos *et al* 2008. Note that `Peacock.jl` doesn't fix the phase of the solutions and your results may differ by a random phase.


### The `Peacock.Zoo` submodule

```julia
using Peacock.Zoo
using Parameters

# Load photonic topological insulator from Wu et al 2015
@unpack geometry, solver, polarisation, G, K, M = make_wu_topo(fourier_space_cutoff)

# Preview the crystal
plot(geometry)
```

> ![Plot of a crystal from the Zoo](figures/example_zoo_geometry.png)

```julia
# Plot the first six bands
figure(figsize=(3,4))
plot_band_diagram(solver, [K,G,M], polarisation,
        bands=1:6, frequency_scale=1/2pi, dk=0.1)
```

> ![Plot of some bands from the Zoo](figures/example_zoo_bands.png)


## Further reading

Photonic crystals - Joannaopoulos
http://ab-initio.mit.edu/book/

Methodology - Rumpf CEM
