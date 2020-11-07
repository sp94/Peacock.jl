module Peacock

using PyPlot
using LinearAlgebra, FFTW

include("utils.jl")

include("geometry.jl")
export Geometry

include("plane_wave_basis.jl")
export BrillouinZoneCoordinate

include("modes.jl")
export Eigenmode, Eigenspace
export overlap

include("solver.jl")
export Solver, TE, TM, solve

include("transforms.jl")
export transform
export symmetry_transform, symmetry_eigvals, symmetry_eigenmodes
export C2, C3, C4, C6, mirror_x, mirror_y

include("plotting.jl")
export plot, get_field

include("band_diagrams.jl")
export plot_band_diagram

include("wilson_loops.jl")
export plot_wilson_loop_winding

# Submodule that defines my commonly used crystals
include("Zoo.jl")

end # module
