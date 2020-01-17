using FFTW

"""
Calculate reciprocal lattice vectors.
"""
function bs(a1, a2)
    b1 = 2pi * [+a2[2],-a2[1]] / (+a1[1]*a2[2]-a1[2]*a2[1])
    b2 = 2pi * [-a1[2],+a1[1]] / (-a2[1]*a1[2]+a2[2]*a1[1])
    @assert dot(a1,b1) ≈ 2pi
    @assert dot(a2,b2) ≈ 2pi
    @assert dot(a1,b2)+1 ≈ 1
    @assert dot(a2,b1)+1 ≈ 1
    return b1, b2
end

"""
Expansion of a geometry onto a basis of spatial harmonics.
"""
struct ConvolvedGeometry
    epc::Array{ComplexF64,2}
    muc::Array{ComplexF64,2}
    a1::Array{Float64,1}
    a2::Array{Float64,1}
    b1::Array{Float64,1}
    b2::Array{Float64,1}
    ps::Array{Int,1}
    qs::Array{Int,1}
    kx::Array{Float64,2}
    ky::Array{Float64,2}
end

"""
Generate convolution matrices.
"""
function convmat(mat, ps::AbstractArray{Int,1}, qs::AbstractArray{Int,1})
    @assert length(ps) == length(qs)
    M = length(ps)
    # Modified from Rumpf CEM Lecture #18
    # "Maxwell's Equations in Fourier Space"
    mat = fftshift(fft(fftshift(mat))) / length(mat)
    p0 = 1 + size(mat,1)÷2
    q0 = 1 + size(mat,2)÷2
    convmat = zeros(ComplexF64, M, M)
    for irow in 1:M, icol in 1:M
        prow, qrow = ps[irow], qs[irow]
        pcol, qcol = ps[icol], qs[icol]
        pfft = prow - pcol
        qfft = qrow - qcol
        convmat[irow,icol] = mat[p0+pfft,q0+qfft]
    end
    # P = length(ps)
    # Q = length(qs)
    # i2s = CartesianIndices((P,Q))
    # for (irow,srow) in enumerate(i2s)
    #     for (icol,scol) in enumerate(i2s)
    #         pfft = ps[srow[1]] - ps[scol[1]]
    #         qfft = qs[srow[2]] - qs[scol[2]]
    #         convmat[irow,icol] = mat[p0+pfft,q0+qfft]
    #     end
    # end
    return convmat
end

function ConvolvedGeometry(g::Geometry, ps::AbstractArray{Int,1}, qs::AbstractArray{Int,1})
    b1, b2 = bs(g.a1, g.a2)
    epc = convmat(g.ep, ps, qs)
    muc = convmat(g.mu, ps, qs)
    # We may as well make kx and ky full matrices now, as we convert later
    # since A_ldiv_B!(A,B) is not defined for sparse B
    kx = diagm(0 => [p*b1[1]+q*b2[1] for (p,q) in zip(ps,qs)])
    ky = diagm(0 => [p*b1[2]+q*b2[2] for (p,q) in zip(ps,qs)])
    return ConvolvedGeometry(epc, muc, g.a1, g.a2, b1, b2, ps, qs, kx, ky)
end

function ConvolvedGeometry(g::Geometry, P::Int, Q::Int)
    @assert isodd(P)
    @assert isodd(Q)
    ps = [p for p in -P÷2:P÷2 for q in -Q÷2:Q÷2]
    qs = [q for p in -P÷2:P÷2 for q in -Q÷2:Q÷2]
    return ConvolvedGeometry(g, ps, qs)
end

function ConvolvedGeometry(g::Geometry, P::Int)
    @assert isodd(P)
    ps, qs = Int[], Int[]
    b1, b2 = bs(g.a1, g.a2)
    @assert norm(b1) ≈ norm(b2) # for now
    d = norm(b1)
    @show d
    for p in -P:P, q in -P:P
        k = p*b1 + q*b2
        if norm(k) <= d*P/2
            push!(ps, p)
            push!(qs, q)
        end
    end
    @show length(ps), length(qs)
    return ConvolvedGeometry(g, ps, qs)
end
