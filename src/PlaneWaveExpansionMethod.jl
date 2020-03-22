module PlaneWaveExpansionMethod

include("brillouin_zones.jl")
export sample_path, BrillouinZoneCoordinate

include("geometry.jl")
export Geometry

include("convolutions.jl")
export ConvolvedGeometry

include("solvers.jl")
export solve, TE, TM

include("postprocessing.jl")
export get_field

include("plotting.jl")
export plot_field, plot_mode

include("brillouin_zones.jl")
export sample_path

include("band_diagrams.jl")
export plot_band_diagram

end # module
