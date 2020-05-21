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
struct Mode
    k0::Vector{Float64}
    data::Vector{ComplexF64}
    weighting::Matrix{ComplexF64}
    basis::PlaneWaveBasis
    label::String
end


"""
Contains eigenvectors in the columns of data
Is always orthonormalised
"""
struct HilbertSpace
    k0::Vector{Float64}
    data::Matrix{ComplexF64}
    weighting::Matrix{ComplexF64}
    basis::PlaneWaveBasis
    function HilbertSpace(k0::Vector{Float64}, data::Matrix{ComplexF64},
                            weighting::Matrix{ComplexF64}, basis::PlaneWaveBasis)
        # Inner constructor guarantees data will be orthonormalised
        data = orthonormalise(data, weighting=weighting)
        return new(k0, data, weighting, basis)
    end
end


"""
To do
"""
function HilbertSpace(modes::Array{Mode,1})
    k0 = modes[1].k0
    data = zeros(ComplexF64, length(modes[1].data), length(modes))
    weighting = modes[1].weighting
    basis = modes[1].basis
    for (col,mode) in enumerate(modes)
        @assert mode.k0 == k0
        @assert mode.weighting == weighting
        @assert mode.basis == basis
        data[:,col] = mode.data
    end
    return HilbertSpace(k0, data, weighting, basis)
end


"""
To do
Also: add test along the lines of:
        hs = HilbertSpace(mode)
        hs = Peacock.shift_k0(hs, 1, 2)
        m = Mode(hs.k0, hs.data[:,1], hs.weighting, hs.basis, "")
        @assert get_field(mode) ≈ get_field(m)
"""
function shift_k0(space::HilbertSpace, dp::Int, dq::Int)
    data = space.data
    ps = space.basis.ps
    qs = space.basis.qs
    data_new = zeros(ComplexF64, size(data))
    for row in 1:size(data,1)
        # Identify the corresponding (p_new,q_new) indices
        # in the plane wave basis
        p_new = ps[row] - dp
        q_new = qs[row] - dq
        row_new = findfirst(x->x==(p_new,q_new), collect(zip(ps,qs)))
        if !isnothing(row_new)
            data_new[row_new,:] = data[row,:]
        end
    end
    k0_new = space.k0 + dp*space.basis.b1 + dq*space.basis.b2
    return HilbertSpace(k0_new, data_new, space.weighting, space.basis)
end
