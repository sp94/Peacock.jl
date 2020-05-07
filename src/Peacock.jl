module Peacock

using PyPlot
using LinearAlgebra, FFTW

include("utils.jl")

include("geometry.jl")
export Geometry

# include("plane_wave_expansion_method/basis.jl")
include("plane_wave_basis.jl")
export Mode, get_field

# include("plane_wave_expansion_method/solver.jl")
include("solver.jl")
export Solver, TE, TM, solve

include("plotting.jl")
export plot

include("band_diagrams.jl")
export BrillouinZoneCoordinate, plot_band_diagram

end # module
