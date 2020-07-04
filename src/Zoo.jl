module Zoo

using Peacock

"""
    wu_ep(x, y, angles, Rs, d1s, d2s, ep_bg, ep_cyl)

Returns the permittivity `ep(x,y)` for photonic crystals similar to Wu et al, 2015.
"""
function wu_ep(x, y, angles, Rs, d1s, d2s, ep_bg, ep_cyl)
    if length(Rs) == 1
        Rs = fill(Rs[1], length(angles))
    end
    if length(d1s) == 1
        d1s = fill(d1s[1], length(angles))
    end
    if length(d2s) == 1
        d2s = fill(d2s[1], length(angles))
    end
    @assert length(Rs) == length(d1s) == length(d2s) == length(angles)
    # Triangular lattice sites
    lx, ly = 1, 1*sqrt(3)
    sites = [
        (0, 0),
        (-lx, 0),
        (+lx, 0),
        (-lx/2, -ly/2),
        (+lx/2, -ly/2),
        (+lx/2, +ly/2),
        (-lx/2, +ly/2),
    ]
    # Hexagonal ring of cylinders at each site
    for (x0,y0) in sites
        for (R,d1,d2,th) in zip(Rs,d1s,d2s,angles)
            x_, y_ = [cosd(th) sind(th); -sind(th) cosd(th)] * [(x-x0); (y-y0)]
            x_ = (x_-R) / (d1/2)
            y_ = (y_-0) / (d2/2)
            if x_^2 + y_^2 <= 1
                return ep_cyl
            end
        end
    end
    return ep_bg
end

"""
    make_wu_generic(Rs, d1s, d2s, polarisation, cutoff;
                angles=[60n for n in 0:5], ep_bg=1, ep_cyl=11.7)

Returns a `NamedTuple` of `geometry`, `solver`, `polarisation`, and the `G`,
`M`, `K` high-symmetry points for photonic crystals similar to Wu et al, 2015.
"""
function make_wu_generic(Rs, d1s, d2s, polarisation, cutoff; angles=[60n for n in 0:5], ep_bg=1, ep_cyl=11.7)
    epf(x,y) = wu_ep(x, y, angles, Rs, d1s, d2s, ep_bg, ep_cyl)
    muf(x,y) = 1.0
    geometry = Geometry(epf, muf, -60, +60, 1/500, 1/500)
    solver = Solver(geometry, cutoff)
    G = BrillouinZoneCoordinate(  0,   0, "Γ")
    M = BrillouinZoneCoordinate(  0, 1/2, "M")
    K = BrillouinZoneCoordinate(1/3, 1/3, "K")
    # Return as a NamedTuple so variables can be easily unpacked using Parameters.jl
    return (; geometry=geometry, solver=solver,
            G=G, M=M, K=K, polarisation=polarisation)
end


"""
    make_wu(a0_div_R, cutoff)

Reproduces "Scheme for Achieving a Topological Photonic Crystal
by Using Dielectric Material", Wu et al, 2015.

Returns a `NamedTuple` of `geometry`, `solver`, `polarisation`, and the
`G`, `M`, `K` high-symmetry points.
"""
function make_wu(a0_div_R, cutoff)
    # Geometry
    R = 1 / a0_div_R # assuming a0 = 1
    d = 2R/3
    return make_wu_generic(R, d, d, TM, cutoff)
end

"""
    make_wu_topo(cutoff::Int)

Reproduces the trivial crystal (a0_div_R=3.125)from  "Scheme for Achieving a
Topological Photonic Crystal by Using Dielectric Material", Wu et al, 2015.

Returns a `NamedTuple` of `geometry`, `solver`, `polarisation`, and the
`G`, `M`, `K` high-symmetry points.
"""
function make_wu_triv(cutoff::Int)
    return make_wu(3.125, cutoff)
end

"""
    make_wu_topo(cutoff::Int)

Reproduces the topological crystal (`a0_div_R=2.9`) from "Scheme for Achieving a
Topological Photonic Crystal by Using Dielectric Material", Wu et al, 2015.

Returns a `NamedTuple` of `geometry`, `solver`, `polarisation`, and the
`G`, `M`, `K` high-symmetry points.
"""
function make_wu_topo(cutoff::Int)
    return make_wu(2.9, cutoff)
end


"""
    make_wu_primitive(cutoff)

Reproduces the primitive unit cell of "Scheme for Achieving a
Topological Photonic Crystal by Using Dielectric Material", Wu et al, 2015

Returns a `NamedTuple` of `geometry`, `solver`, `polarisation`, and the
`G`, `M`, `K` high-symmetry points.
"""
function make_wu_primitive(cutoff)
    # Geometry
    R = 1 / 3 # assuming a0 = 1
    d = 2R/3
    x = 1 / sqrt(3)
    a1 = [x*cosd(-30), x*sind(-30)]
    a2 = [x*cosd(+30), x*sind(+30)]
    epf(x,y) = wu_ep(x, y, R, d, d)
    muf(x,y) = 1.0
    g = Geometry(epf, muf, a1, a2, 1/500, 1/500)
    solver = Solver(cg, TM)
    # BZ symmetry points
    G = BrillouinZoneCoordinate(  0,   0, "Γ")
    M = BrillouinZoneCoordinate(  0, 1/2, "M")
    K = BrillouinZoneCoordinate(1/3, 1/3, "K")
    return g, solver, G, M, K
end


"""
    make_dePaz(d1::Real, d2::Real, cutoff::Int; R::Real=1/3)

Reproduces "Engineering fragile topology in photonic crystals:
Topological quantum chemistry of light", Blanco de Paz et al, 2019

Returns a `NamedTuple` of `geometry`, `solver`, `polarisation`, and the
`G`, `M`, `K` high-symmetry points.
"""
function make_dePaz(d1::Real, d2::Real, cutoff::Int; R::Real=1/3)
    return make_wu_generic(R, d1, d2, TM, cutoff)
end

"""
    make_dePaz_triv(cutoff::Int)

Reproduces the trivial (`d1=0.52`, `d2=0.31`) crystal from "Engineering fragile
topology in photonic crystals: Topological quantum chemistry of light",
Blanco de Paz et al, 2019

Returns a `NamedTuple` of `geometry`, `solver`, `polarisation`, and the
`G`, `M`, `K` high-symmetry points.
"""
function make_dePaz_triv(cutoff::Int)
    return make_dePaz(0.52, 0.31, cutoff)
end

"""
    make_dePaz_frag(cutoff::Int)

Reproduces the trivial (`d1=0.4`, `d2=0.13`) crystal from "Engineering fragile
topology in photonic crystals: Topological quantum chemistry of light",
Blanco de Paz et al, 2019

Returns a `NamedTuple` of `geometry`, `solver`, `polarisation`, and the
`G`, `M`, `K` high-symmetry points.
"""
function make_dePaz_frag(cutoff::Int)
    return make_dePaz(0.4, 0.13, cutoff)
end

"""
    make_dePaz_OAL(cutoff::Int)

Reproduces the obstructed atomic limit (`d1=0.4`, `d2=0.61`) crystal from
"Engineering fragile topology in photonic crystals: Topological quantum
chemistry of light", Blanco de Paz et al, 2019

Returns a `NamedTuple` of `geometry`, `solver`, `polarisation`, and the
`G`, `M`, `K` high-symmetry points.
"""
function make_dePaz_OAL(cutoff::Int)
    return make_dePaz(0.4, 0.61, cutoff)
end


"""
    make_zhu(r, l, cutoff; angles=[60n for n in 0:5])

Reproduces "Topological transitions in continuously deformed photonic crystals"
Zhu et al, 2018.

- `r` is the radius of the cylindrical air holes.
- `l` is the radius of the hexagonal ring that the air holes form

Note that 'l' in Zhu et al is the same as 'R' in wu_ep.

Returns a `NamedTuple` of `geometry`, `solver`, `polarisation`, and the
`G`, `M`, `K` high-symmetry points.
"""
function make_zhu(r, l, cutoff; angles=[60n for n in 0:5])
    return make_wu_generic(l, 2r, 2r, TE, cutoff, ep_cyl=1, ep_bg=10.2)
end


"""
    make_xu(r1, r2, cutoff; wedges=[(0,360)])

Reproduces "Accidental degeneracy in photonic bands and topological phase
transitions in two-dimensional core-shell dielectric photonic crystals",
Xu et al, 2016.

The rings can be split using `wedges = [(start,stop), (start,stop)...]`.

Returns a `NamedTuple` of `geometry`, `solver`, `polarisation`, and the
`G`, `M`, `K` high-symmetry points.
"""
function make_xu(r1, r2, cutoff; wedges=[(0,360)])
    function epf(x,y)
        lx, ly = 1, 1*sqrt(3)
        sites = [
            (    0,     0),
            (-lx/2, -ly/2),
            (+lx/2, -ly/2),
            (+lx/2, +ly/2),
            (-lx/2, +ly/2),
        ]
        for (x0,y0) in sites
            x_, y_ = x-x0, y-y0
            if r1^2 <= x_^2 + y_^2 <= r2^2
                θ = atan(y_,x_) * 180/pi
                for (θ1,θ2) in wedges
                    dθ = mod(θ - θ1, 360)
                    if dθ <= (θ2-θ1)
                        return 12.0
                    end
                end
            end
        end
        return 1.0
    end
    muf(x,y) = 1.0
    geometry = Geometry(epf, muf, -60, +60, 1/500, 1/500)
    solver = Solver(geometry, cutoff)
    G = BrillouinZoneCoordinate(  0,   0, "Γ")
    M = BrillouinZoneCoordinate(  0, 1/2, "M")
    K = BrillouinZoneCoordinate(1/3, 1/3, "K")
    return (; geometry=geometry, solver=solver,
            G=G, M=M, K=K, polarisation=TM)
end


export make_wu, make_wu_triv, make_wu_topo, make_wu_primitive
export make_dePaz, make_dePaz_triv, make_dePaz_frag, make_dePaz_OAL
export make_zhu
export make_xu

end # module
