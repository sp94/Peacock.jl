function round_with_tolerance(X, tol=1e-6)
    X_ = Int.(round.(X))
    @assert all(abs.(X-X_) .< 1e-6)
    return X_
end

function transform(data::AbstractArray{ComplexF64}, basis::PlaneWaveBasis, k_map::Function)
    # Identify how the basis transforms
    b1 = basis.b1
    b2 = basis.b2
    B = [b1 b2]
    shift = k_map([0,0])
    dp, dq = round_with_tolerance(B \ shift)
    b1_ = k_map(b1) - shift
    b2_ = k_map(b2) - shift
    X = round_with_tolerance(B \ [b1_ b2_])
    p_map(p,q) = X[1,1]*p + X[1,2]*q - dp
    q_map(p,q) = X[2,1]*p + X[2,2]*q - dq
    # Create transformed copy of the data
    data_new = zeros(ComplexF64, size(data))
    for row in 1:size(data,1)
        p = basis.ps[row]
        q = basis.qs[row]
        p_ = p_map(p,q)
        q_ = q_map(p,q)
        row_ = findfirst(x->x==(p_,q_), collect(zip(basis.ps,basis.qs)))
        if !isnothing(row_)
            data_new[row_,:] = data[row,:]
        end
    end
    return data_new
end

function transform(mode::Eigenmode, k_map::Function)
    k0_ = k_map(mode.k0)
    data_ = transform(mode.data, mode.basis, k_map)
    return Eigenmode(k0_, mode.frequency, data_, mode.weighting, mode.basis, mode.label)
end

function transform(space::Eigenspace, k_map::Function)
    k0_ = k_map(space.k0)
    data_ = transform(space.data, space.basis, k_map)
    Eigenspace(k0_, data_, space.weighting, space.basis)
end


function eigen(space::Eigenspace, k_map::Function)
    # shift such that _k_map(space.k0) ≈ space.k0
    _k_map(k) = k_map(k) + space.k0 - k_map(space.k0)
    space_ = transform(space, _k_map)
    xlaps = overlaps(space, space_)
    vals, vecs = eigen(xlaps)
    modes = Eigenmode[]
    for n in 1:length(vals)
        data = space.data * vecs[:,n]
        mode = Eigenmode(space.k0, vals[n], data, space.weighting, space.basis, "")
        push!(modes, mode)
    end
    return modes
end

function eigen(mode::Eigenmode, k_map::Function)
    space = Eigenspace(mode)
    return eigen(space, k_map)
end

# Useful transformations
C2(k) = rotation_matrix(180)*k
C3(k) = rotation_matrix(120)*k
C4(k) = rotation_matrix( 90)*k
C6(k) = rotation_matrix( 60)*k
mirror_x(k) = [-k[1], +k[2]]
mirror_y(k) = [+k[1], -k[2]]
