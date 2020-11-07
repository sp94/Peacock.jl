"""
    Eigenmode(k0, frequency, data, weighting, basis, label)

Eigenmode of a photonic crystal expressed on a plane-wave `basis` with a weighted inner product.
"""
struct Eigenmode
    k0::Vector{Float64}
    data::Vector{ComplexF64}
    weighting::Matrix{ComplexF64}
    basis::PlaneWaveBasis
    eigenvalue::ComplexF64
    data_label::String
    eigenvalue_label::String
end

function overlap(m1::Eigenmode, m2::Eigenmode)
    LHS = m1.data
    RHS = m2.weighting * m2.data
    return dot(LHS, RHS)
end


"""
    Eigenspace(k0, data, weighting, basis)

eigenspace spanned by the eigenvectors in each column of `data`.

The eigenvectors will be orthonormalised using Gram-Schmidt orthonormalisation,
see [`orthonormalise`](@ref).
"""
struct Eigenspace
    k0::Vector{Float64}
    data::Matrix{ComplexF64}
    weighting::Matrix{ComplexF64}
    basis::PlaneWaveBasis
    data_label::String
    function Eigenspace(k0::Vector{Float64}, data::Matrix{ComplexF64},
                            weighting::Matrix{ComplexF64}, basis::PlaneWaveBasis,
                            data_label::String="")
        # Inner constructor guarantees data will be orthonormalised
        data = orthonormalise(data, weighting=weighting)
        return new(k0, data, weighting, basis, data_label)
    end
end


"""
    Eigenspace(mode::Eigenmode)

Returns the eigenspace spanned by `mode`.

The `data` of the eigenspace is guaranteed to be normalised
under the weighting of the `mode`.
"""
function Eigenspace(mode::Eigenmode)
    return Eigenspace([mode])
end



"""
    Eigenspace(modes::Array{Eigenmode,1})

Returns the eigenspace spanned by the `modes`.

The `data` of the eigenspace is guaranteed to be orthonormal
under the weighting of the `modes`.
"""
function Eigenspace(modes::Array{Eigenmode,1})
    k0 = modes[1].k0
    data = zeros(ComplexF64, length(modes[1].data), length(modes))
    weighting = modes[1].weighting
    basis = modes[1].basis
    data_label = modes[1].data_label
    for (col,mode) in enumerate(modes)
        @assert mode.k0 == k0
        @assert mode.weighting == weighting
        @assert mode.basis == basis
        @assert mode.data_label == data_label
        data[:,col] = mode.data
    end
    return Eigenspace(k0, data, weighting, basis, data_label)
end


"""
    normalise(data; weighting=I)

Normalisation of vectors with a weighted inner product.
"""
function normalise(data; weighting=I)
    out = copy(data)
    for n in 1:size(out,2)
        out[:,n] /= sqrt(abs(dot(out[:,n],weighting*out[:,n])))
    end
    return out
end


"""
    orthonormalise(data; weighting=I)

Gram-Schmidt orthonormalisation of vectors with a weighted inner product.
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
