"""
    Geometry(a1::Array{Real,1}, b1::Array{Real,1}, ep::Array{ComplexF64,2}, mu::Array{ComplexF64,2})

2D geometry defined in real space, with lattice vectors `a1` and `a2`,
and relative permeability and permittivity `ep` and `mu`, respectively.
"""
struct Geometry
    a1::Array{Real,1}
    a2::Array{Real,1}
    ep::Array{ComplexF64,2}
    mu::Array{ComplexF64,2}
end


"""
    Geometry(epf::Function, muf::Function, a1::Array{<:Real,1}, a2::Array{<:Real,1}, d1::Real, d2::Real)

Generate geometry with permittivity `epf(x,y)` and permeability `muf(x,y)`.

The real space lattice vectors, `a1` and `a2`, define the unit cell.
The grid resolution along each lattice vector is `d1` and `d2`, respectively.
"""
function Geometry(epf::Function, muf::Function,
            a1::Array{<:Real,1}, a2::Array{<:Real,1}, d1::Real, d2::Real)
    P = ceil(Int, norm(a1)/d1)
    Q = ceil(Int, norm(a2)/d2)
    ps = range(-0.5, stop=0.5, length=P)[1:end-1]
    qs = range(-0.5, stop=0.5, length=Q)[1:end-1]
    ps = ps .+ step(ps)/2 # ensure that p and q are centered on zero
    qs = qs .+ step(qs)/2
    xs = [p*a1[1]+q*a2[1] for p in ps, q in qs]
    ys = [p*a1[2]+q*a2[2] for p in ps, q in qs]
    ep = epf.(xs, ys)
    mu = muf.(xs, ys)
    return Geometry(a1, a2, ep, mu)
end


"""
    Geometry(epf::Function, muf::Function, a1::Array{<:Real,1}, a2::Array{<:Real,1}, d1::Real, d2::Real)

Generate geometry with permittivity `epf(x,y)` and permeability `muf(x,y)`.

The real space lattice vectors are assumed to have unit length and are at
angles `a1_deg` and `a2_deg` counter-clockwise from the x-axis.
The grid resolution along each lattice vector is `d1` and `d2`, respectively.
"""
function Geometry(epf::Function, muf::Function, a1_deg::Real, a2_deg::Real, d1::Real, d2::Real)
    a1 = [cosd(a1_deg), sind(a1_deg)]
    a2 = [cosd(a2_deg), sind(a2_deg)]
    return Geometry(epf, muf, a1, a2, d1, d2)
end
