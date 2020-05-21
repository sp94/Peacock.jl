"""
To do
"""
struct PlaneWaveBasis
    b1::Vector{Float64}
    b2::Vector{Float64}
    ps::Vector{Int}
    qs::Vector{Int}
    kxs::Vector{ComplexF64}
    kys::Vector{ComplexF64}
end


"""
To do
Useful for shifting basis or for including a Bloch phase
"""
function PlaneWaveBasis(basis::PlaneWaveBasis, dk::AbstractVector{<:Real})
    kxs = [kx+dk[1] for kx in basis.kxs]
    kys = [ky+dk[2] for ky in basis.kys]
    return PlaneWaveBasis(basis.b1, basis.b2, basis.ps, basis.qs, kxs, kys)
end


"""
To do
"""
function PlaneWaveBasis(geometry::Geometry, cutoff::Int)
    @assert isodd(cutoff)
    b1, b2 = as_to_bs(geometry.a1, geometry.a2)
    @assert norm(b1) ≈ norm(b2) # for now
    ps, qs = Int[], Int[]
    cutoff_radius = norm(b1) * cutoff / 2
    for p in -cutoff:cutoff, q in -cutoff:cutoff
        k = p*b1 + q*b2
        if norm(k) <= cutoff_radius
            push!(ps, p)
            push!(qs, q)
        end
    end
    kxs = [p*b1[1]+q*b2[1] for (p,q) in zip(ps,qs)]
    kys = [p*b1[2]+q*b2[2] for (p,q) in zip(ps,qs)]
    return PlaneWaveBasis(b1, b2, ps, qs, kxs, kys)
end


"""
To do
"""
struct BrillouinZoneCoordinate
    p::Float64
    q::Float64
    label::String
end

function BrillouinZoneCoordinate(p::Real, q::Real; label="")
    return BrillouinZoneCoordinate(p, q, "")
end


"""
To do
"""
function get_k(c::BrillouinZoneCoordinate, basis::PlaneWaveBasis)
    return c.p*basis.b1 + c.q*basis.b2
end


"""
To do
"""
function get_field(us::AbstractVector{<:Complex}, basis::PlaneWaveBasis;
                    k0=[0,0], t1s=-0.5:0.01:0.5, t2s=-0.5:0.01:0.5)
    # Generate real space coordinates
    b1, b2 = basis.b1, basis.b2
    a1, a2 = bs_to_as(b1, b2)
    xs = [t1*a1[1]+t2*a2[1] for t1 in t1s, t2 in t2s]
    ys = [t1*a1[2]+t2*a2[2] for t1 in t1s, t2 in t2s]
    out = zeros(ComplexF64, length(t1s), length(t2s))
    # Inverse Fourier transform
    kxs = k0[1] .+ basis.kxs
    kys = k0[2] .+ basis.kys
    for (kx,ky,u) in zip(kxs, kys, us)
        for c in CartesianIndices(out)
            phase = kx*xs[c] + ky*ys[c]
            out[c] += u * exp(1im*phase)
        end
    end
    return out
end


"""
To do
"""
struct Mode <: AbstractVector{ComplexF64}
    k0::Vector{Float64}
    data::Vector{ComplexF64}
    weighting::Matrix{ComplexF64}
    basis::PlaneWaveBasis
    label::String
end
Base.size(A::Mode) = size(A.data)
Base.getindex(A::Mode, I::Int) = A.data[I]
