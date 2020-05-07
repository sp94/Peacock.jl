"""
To do
"""
@enum Polarisation TE TM


"""
To do
"""
struct Solver
    basis::PlaneWaveBasis
    epc::Matrix{ComplexF64}
    muc::Matrix{ComplexF64}
end


"""
To do
"""
function Solver(geometry::Geometry, cutoff::Int)
    basis = PlaneWaveBasis(geometry, cutoff)
    epc = convmat(geometry.ep, basis)
    muc = convmat(geometry.mu, basis)
    return Solver(basis, epc, muc)
end


"""
Generate convolution matrices, see Raymond Rumpf's CEM Lecture #18,
"Maxwell's Equations in Fourier Space", for further reading.
"""
function convmat(mat, basis::PlaneWaveBasis)
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
To do
"""
function solve(solver::Solver, k::AbstractVector{<:Real}, polarisation::Polarisation; bands=:)
    # Get left and right hand sides of the generalised eigenvalue problem
    basis = solver.basis
    bloch_basis = PlaneWaveBasis(basis, k)
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
    modes = [Mode(k,data,weighting,basis,label) for data in eachcol(modes_data)]
    # modes = orthonormalise(modes)
    return freqs, modes
end


"""
To do
"""
function solve(solver::Solver, x::BrillouinZoneCoordinate, polarisation::Polarisation, bands=:)
    k = get_k(x, solver.basis)
    return solve(solver, k, polarisation, bands=bands)
end
