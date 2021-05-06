using PyPlot
using Peacock

function epf(x,y)
    # equation of a circle with radius 0.2a
    if x^2+y^2 <= 0.2^2
        # dielectric inside the circle
        return 8.9
    else
        # air outside the circle
        return 1
    end
end

# Permeability is unity everywhere
function muf(x,y)
    return 1
end

a1 = [1, 0]  # first lattice vector
a2 = [0, 1]  # second lattice vector
d1 = 0.01  # resolution along first lattice vector
d2 = 0.01  # resolution along second lattice vector
geometry = Geometry(epf, muf, a1, a2, d1, d2)
fourier_space_cutoff = 9
solver_CPU = Solver(geometry, fourier_space_cutoff)
solver_GPU = Solver(geometry, fourier_space_cutoff, GPU=true)
G = BrillouinZoneCoordinate(  0,   0, "Γ")
X = BrillouinZoneCoordinate(1/2,   0, "X")
M = BrillouinZoneCoordinate(1/2, 1/2, "M")
ks = [G,X,M,G]

@show fourier_space_cutoff
@show typeof(solver_CPU.epc)
@show typeof(solver_GPU.muc)


function test_bands(solver, dk)
    figure(figsize=(4,3))
    plot_band_diagram(solver, ks, TE, color="red",
                bands=1:4, dk=dk, frequency_scale=1/2pi)
    plot_band_diagram(solver, ks, TM, color="blue",
                bands=1:4, dk=dk, frequency_scale=1/2pi)
    ylim(0,0.8)
end

# call once to make sure functions are compiled
test_bands(solver_CPU, 2)
test_bands(solver_GPU, 2)

close("all")
@time test_bands(solver_CPU, 0.1)
@time test_bands(solver_GPU, 0.1)
show()
