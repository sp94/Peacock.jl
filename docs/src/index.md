![](assets/logo.png)

# Introduction

`Peacock.jl` - or the **P**lane-wave **E**xpansion **A**pproach to **C**haracterising **O**ptical **C**rystals in **k**-space - is a Julia package for studying photonic crystals using the Plane Wave Expansion Method.


* Solve for...
  * Transverse electric (TE) and transverse magnetic (TM) modes of 2D photonic crystals
  * Non-orthogonal unit cells
  * Inhomogeneous permittivity and/or permeability
  * Chern numbers of topological photonic crystals using built-in Wilson loop methods


* Focused on ease of use
  * Install with one line in Julia's package manager
  * Simple visualisation of geometry, fields, and fully labelled band diagrams


Photonic crystals are materials whose optical properties arise from the structuring of the material when the size of the structures are comparable to the wavelengths of light. `Peacock.jl` is named for the irridescent colours of peacock feathers which arise not from pigmentation but from their photonic crystal structure, as shown below.

> ![Different zooms of a Peacock](assets/peacock_feathers_zoom.png)
> Peacock by allanlau2000 from pixabay. Feather by suju from pixabay. Electron microscope image of photonic crystal structure from Zi, Jian, et al. "Coloration strategies in peacock feathers.",  *Proceedings of the National Academy of Sciences* 100.22 (2003): 12576-12578. Copyright (2003) National Academy of Sciences.

As well as occuring naturally as in animals such as peacocks, advances in nanofabrication mean that 'designer' photonic crystals can now be manufactured for unprecedented control over the flow of light, with applications ranging from optical fibers to photonic circuitry. Photonic crystals are also a promising platform for more exotic materials like topological insulators.


## Index
```@contents
Pages = [
        "tutorials/getting_started.md",
        "reference/index.md"]
```
