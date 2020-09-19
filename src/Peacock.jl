module Peacock

using PyPlot
using LinearAlgebra, FFTW

include("utils.jl")

include("geometry.jl")
export Geometry

include("plane_wave_basis.jl")
export BrillouinZoneCoordinate

include("modes.jl")
export Mode, HilbertSpace, get_field

include("solver.jl")
export Solver, TE, TM, solve

include("transforms.jl")
export transform

include("plotting.jl")
export plot

include("band_diagrams.jl")
export plot_band_diagram

include("wilson_loops.jl")
export plot_wilson_loop_winding

# Submodule that defines my commonly used crystals
include("Zoo.jl")

end # module
