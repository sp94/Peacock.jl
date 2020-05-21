"""
Calculate reciprocal lattice vectors from real space lattice vectors.
"""
function as_to_bs(a1, a2)
    b1 = 2pi * [+a2[2],-a2[1]] / (+a1[1]*a2[2]-a1[2]*a2[1])
    b2 = 2pi * [-a1[2],+a1[1]] / (-a2[1]*a1[2]+a2[2]*a1[1])
    @assert dot(a1,b1) ≈ 2pi
    @assert dot(a2,b2) ≈ 2pi
    @assert dot(a1,b2)+1 ≈ 1
    @assert dot(a2,b1)+1 ≈ 1 # TO DO - move assertion into test
    return b1, b2
end


"""
Calculate real space lattice vectors from reciprocal lattice vectors.
"""
function bs_to_as(b1, b2)
    # Actually the same formula, but I think this naming makes the code clearest
    # TO DO - implement as test: bs_to_as(as_to_bs(a1,a2)) == a1,a2
    return as_to_bs(b1, b2)
end


"""
A sparse diagonal matrix that can be used in left division (D \\ X)
"""
struct DiagonalMatrix <: AbstractMatrix{ComplexF64}
    diag::AbstractVector{ComplexF64}
end
Base.size(A::DiagonalMatrix) = (length(A.diag), length(A.diag))
Base.getindex(A::DiagonalMatrix, I::Vararg{Int,2}) = I[1]==I[2] ? A.diag[I[1]] : 0


"""
To do
"""
function normalise(data; weighting=I)
    out = copy(data)
    for n in 1:size(out,2)
        out[:,n] /= sqrt(abs(dot(out[:,n],weighting*out[:,n])))
    end
    return out
end


"""
To do
"""
function orthonormalise(data; weighting=I)
    out = copy(data)
    out[:,1] = normalise(out[:,1], weighting=weighting)
    # Gram-Schmidt orthogonalisation
    projector = zeros(size(out,1),size(out,1))
    for n in 2:size(out,2)
        projector += out[:,n-1]*out[:,n-1]' * weighting
        out[:,n] = out[:,n] - projector*out[:,n]
        out[:,n] = normalise(out[:,n], weighting=weighting)
    end
    return out
end
