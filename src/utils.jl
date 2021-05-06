"""
    as_to_bs(a1, a2)

Calculate reciprocal lattice vectors, `b1` and `b2`, from the real space lattice
vectors, `a1` and `a2`.
"""
function as_to_bs(a1, a2)
    b1 = 2pi * [+a2[2],-a2[1]] / (+a1[1]*a2[2]-a1[2]*a2[1])
    b2 = 2pi * [-a1[2],+a1[1]] / (-a2[1]*a1[2]+a2[2]*a1[1])
    return b1, b2
end


"""
    bs_to_as(b1, b2)

Calculate real space lattice vectors, `a1` and `a2`, from reciprocal lattice
vectors, `b1` and `b2`.

This is actually the same as `as_to_bs`, but I think having both functions makes
the intention of the code more obvious.
"""
function bs_to_as(b1, b2)
    return as_to_bs(b1, b2)
end


# """
#     DiagonalMatrix(diag::AbstractVector{ComplexF64})

# A sparse diagonal matrix that can be used in left division (D \\ X)
# """
# struct DiagonalMatrix <: AbstractMatrix{ComplexF64}
#     diag::AbstractVector{ComplexF64}
# end
# Base.size(A::DiagonalMatrix) = (length(A.diag), length(A.diag))
# Base.getindex(A::DiagonalMatrix, I::Vararg{Int,2}) = I[1]==I[2] ? A.diag[I[1]] : 0
