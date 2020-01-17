"""
2D geometry defined in real space.
"""
struct Geometry
    epf::Function
    muf::Function
    a1::Array{Real,1}
    a2::Array{Real,1}
    ep::Array{ComplexF64,2}
    mu::Array{ComplexF64,2}
end

"""
Generate geometry on real-space grid from the permittivity and permeability
functions, epf(x,y) and muf(x,y), respectively.

The real space lattice vectors, a1 and a2, define the unit cell.
The grid resolution along each lattice vector is d1 and d2, respectively.
"""
function Geometry(epf::Function, muf::Function, a1::Array{<:Real,1}, a2::Array{<:Real,1}, d1::Real, d2::Real)
    ps = range(-0.5, stop=0.5, length=ceil(Int,norm(a1)/d1))[1:end-1]
    qs = range(-0.5, stop=0.5, length=ceil(Int,norm(a2)/d2))[1:end-1]
    ps = ps .+ step(ps)/2 # ensure that p and q are centered on zero
    qs = qs .+ step(qs)/2
    xs = [p*a1[1]+q*a2[1] for p in ps, q in qs]
    ys = [p*a1[2]+q*a2[2] for p in ps, q in qs]
    ep = epf.(xs, ys)
    mu = muf.(xs, ys)
    return Geometry(epf, muf, a1, a2, ep, mu)
end

function Geometry(epf::Function, muf::Function, a1_deg::Real, a2_deg::Real, d1::Real, d2::Real)
    a1 = [cosd(a1_deg), sind(a1_deg)]
    a2 = [cosd(a2_deg), sind(a2_deg)]
    return Geometry(epf, muf, a1, a2, d1, d2)
end
