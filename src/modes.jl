"""
    Mode(k0, frequency, data, weighting, basis, label)

Eigenmode of a photonic crystal expressed on a plane-wave `basis` with a weighted inner product.
"""
struct Mode
    k0::Vector{Float64}
    frequency::ComplexF64
    data::Vector{ComplexF64}
    weighting::Matrix{ComplexF64}
    basis::PlaneWaveBasis
    label::String
end


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
    HilbertSpace(modes::Array{Mode,1})

Returns the Hilbert space spanned by the `modes`.

The `data` of the Hilbert space is guaranteed to be orthonormal
under the weighting of the `modes`.
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
