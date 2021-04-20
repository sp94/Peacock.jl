"""
Transverse electric (TE) or transverse magnetic (TM) polarisation.
"""
@enum Polarisation TE TM

"""
    Solver(basis::PlaneWaveBasis, epc::Matrix{ComplexF64}, muc::Matrix{ComplexF64})

A plane-wave expansion method solver, where `epc` and `muc` are the convolution
matrices of the permittivity and permeability, respectively, for the given
`basis` of plane waves.
"""
struct Solver
    basis::PlaneWaveBasis
    epc::Matrix{ComplexF64}
    muc::Matrix{ComplexF64}
end


"""
    Solver(geometry::Geometry, basis::PlaneWaveBasis)

Approximate the geometry using the given `basis` of plane waves.
"""
function Solver(geometry::Geometry, basis::PlaneWaveBasis)
    epc = convmat(geometry.ep, basis)
    muc = convmat(geometry.mu, basis)
    return Solver(basis, epc, muc)
end


"""
    Solver(geometry::Geometry, cutoff::Int)

Approximate the geometry using a basis of plane waves truncated in a circle.

The circle has a diameter of `cutoff` Brillouin zones. Increasing the `cutoff`
will increase the number of plane waves leading to a more accurate solution.
It is assumed that `norm(b1) == norm(b2)`.
"""
function Solver(geometry::Geometry, cutoff::Int)
    basis = PlaneWaveBasis(geometry, cutoff)
    return Solver(geometry, basis)
end


"""
    Solver(geometry::Geometry, cutoff_b1::Int, cutoff_b2::Int)

Approximate the geometry using a basis of plane waves truncated in a rhombus.

The rhombus has lengths `cutoff_b1` and `cutoff_b2` in the `b1` and `b2`
directions, respectively.
"""
function Solver(geometry::Geometry, cutoff_b1::Int, cutoff_b2::Int)
    basis = PlaneWaveBasis(geometry, cutoff_b1, cutoff_b2)
    return Solver(geometry, basis)
end


"""
    convmat(mat::AbstractMatrix, basis::PlaneWaveBasis)

Generate convolution matrices, see Raymond Rumpf's CEM Lecture #18,
"Maxwell's Equations in Fourier Space", for further reading.
"""
function convmat(mat::AbstractMatrix, basis::PlaneWaveBasis)
    @assert length(basis.ps) == length(basis.qs)
    M = length(basis.ps)
    # Modified from Rumpf CEM Lecture #18
    # "Maxwell's Equations in Fourier Space"
    mat = fftshift(fft(fftshift(mat))) / length(mat)
    p0 = 1 + size(mat,1)÷2
    q0 = 1 + size(mat,2)÷2
    convmat = zeros(ComplexF64, M, M)
    for irow in 1:M, icol in 1:M
        prow, qrow = basis.ps[irow], basis.qs[irow]
        pcol, qcol = basis.ps[icol], basis.qs[icol]
        pfft = prow - pcol
        qfft = qrow - qcol
        convmat[irow,icol] = mat[p0+pfft,q0+qfft]
    end
    return convmat
end


"""
    solve(solver::Solver, k::AbstractVector{<:Real}, polarisation::Polarisation, GPU::Int; bands=:)

Calculate the eigenmodes of a photonic crystal at position `k` in reciprocal space. Add GPU support
"""
function solve(solver::Solver, k::AbstractVector{<:Real}, polarisation::Polarisation, GPU::Int; bands=:)
    # Get left and right hand sides of the generalised eigenvalue problem
    basis = solver.basis
    if GPU == 1
        epc, muc = CuArray(solver.epc), CuArray(solver.muc)
        Kx = CuArray(DiagonalMatrix(basis.kxs) + k[1]*I)
        Ky = CuArray(DiagonalMatrix(basis.kys) + k[2]*I)
        if polarisation == TE
            LHS = Kx/epc*Kx + Ky/epc*Ky
            RHS = muc
            label = L"H_z"
        elseif polarisation == TM
            LHS = Kx/muc*Kx + Ky/muc*Ky
            RHS = epc
            label = L"E_z"
        end
        # Sometimes the generalised eigenvalue problem solver
        # fails near Γ when the crystals are symmetric.
        # In these cases, rewrite as a regular eigenvalue problem
        freqs_squared, modes_data = eigen(Array(RHS \ LHS));
    else
        epc, muc = solver.epc, solver.muc
        Kx = DiagonalMatrix(basis.kxs) + k[1]*I
        Ky = DiagonalMatrix(basis.kys) + k[2]*I
        if polarisation == TE
            LHS = Kx/epc*Kx + Ky/epc*Ky
            RHS = muc
            label = L"H_z"
        elseif polarisation == TM
            LHS = Kx/muc*Kx + Ky/muc*Ky
            RHS = epc
            label = L"E_z"
        end
        # Sometimes the generalised eigenvalue problem solver
        # fails near Γ when the crystals are symmetric.
        # In these cases, rewrite as a regular eigenvalue problem
        freqs_squared, modes_data = try
            eigen(LHS, RHS)
        catch
            eigen(RHS \ LHS)
        end
    end
    
    freqs = sqrt.(freqs_squared)
    # Sort by increasing frequency
    idx = sortperm(freqs, by=real)
    freqs = freqs[idx][bands]
    modes_data = modes_data[:,idx][:,bands]
    # Eigenmodes are weighted by the RHS of the generalised eigenvalue problem
    weighting = RHS
    modes = Mode[]
    for i in 1:length(freqs)
        mode = Mode(k, freqs[i], modes_data[:,i], weighting, basis, label)
        push!(modes, mode)
    end
    return modes
end

function solve(solver::Solver, k::AbstractVector{<:Real}, polarisation::Polarisation; bands=:)
    # Get left and right hand sides of the generalised eigenvalue problem
    basis = solver.basis
    epc, muc = solver.epc, solver.muc
    Kx = DiagonalMatrix(basis.kxs) + k[1]*I
    Ky = DiagonalMatrix(basis.kys) + k[2]*I
    if polarisation == TE
        LHS = Kx/epc*Kx + Ky/epc*Ky
        RHS = muc
        label = L"H_z"
    elseif polarisation == TM
        LHS = Kx/muc*Kx + Ky/muc*Ky
        RHS = epc
        label = L"E_z"
    end
    # Sometimes the generalised eigenvalue problem solver
    # fails near Γ when the crystals are symmetric.
    # In these cases, rewrite as a regular eigenvalue problem
    freqs_squared, modes_data = try
        eigen(LHS, RHS)
    catch
        eigen(RHS \ LHS)
    end
    
    freqs = sqrt.(freqs_squared)
    # Sort by increasing frequency
    idx = sortperm(freqs, by=real)
    freqs = freqs[idx][bands]
    modes_data = modes_data[:,idx][:,bands]
    # Eigenmodes are weighted by the RHS of the generalised eigenvalue problem
    weighting = RHS
    modes = Mode[]
    for i in 1:length(freqs)
        mode = Mode(k, freqs[i], modes_data[:,i], weighting, basis, label)
        push!(modes, mode)
    end
    return modes
end


"""
    solve(solver::Solver, x::BrillouinZoneCoordinate, polarisation::Polarisation; bands=:)

Calculate the eigenmodes of a photonic crystal at position `k=x.p*b1 + x.q*b2` in reciprocal space.
"""
function solve(solver::Solver, x::BrillouinZoneCoordinate, polarisation::Polarisation, bands=:)
    k = get_k(x, solver.basis)
    return solve(solver, k, polarisation, bands=bands)
end
