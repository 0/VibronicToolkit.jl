# Diagonal density matrix.

"""
Density matrix that is diagonal in position and off-diagonal in surfaces.
"""
abstract type AbstractDiagonalDensity end

"""
Diagonal density for a small system at finite temperature.
"""
struct DiagonalDensity <: AbstractDiagonalDensity
    "Diagonal density."
    density::Array{Float64,4}
end

"""
    DiagonalDensity(sys::System{S,2}, beta::Float64, basis_size::Int, extent::NTuple{4,Float64}; lengths::NTuple{2,Int}=(101, 101), progress_output::IO=stderr)

Compute the diagonal (in position) density matrix for the 2-mode system `sys`
at `beta` using `basis_size` basis functions.

The extent for each mode is given in `extent`, and the number of points is
given in `lengths`.

The progress meter is written to `progress_output`.
"""
function DiagonalDensity(sys::System{S,2}, beta::Float64, basis_size::Int, extent::NTuple{4,Float64}; lengths::NTuple{2,Int}=(101, 101), progress_output::IO=stderr) where {S}
    basis = Basis(sys, basis_size)
    h0, V = operators(basis, sys)
    basis_vectors = vectors(basis)

    F = eigen(Symmetric(h0 + V))
    Es = F.values
    Vs = F.vectors

    Z = sum(exp.(-beta * Es))

    M = 2
    q1_min, q1_max, q2_min, q2_max = extent
    q1_length, q2_length = lengths
    density = zeros(q2_length, q1_length, S, S)
    wfs = Array{Float64}(undef, basis.dim, S)

    meter = Progress(q1_length, output=progress_output)
    for (i, q1) in ProgressWrapper(enumerate(range(q1_min; stop=q1_max, length=q1_length)), meter)
        for (j, q2) in enumerate(range(q2_min; stop=q2_max, length=q2_length))
            fill!(wfs, 0.0)
            for idx in 1:basis.dim
                s = basis_vectors[M+1, idx]
                wfs[:, s] .+= ho_wf(basis_vectors[1:M, idx], [q1, q2]) .* Vs[idx, :]
            end
            for s1 in 1:S
                for s2 in 1:S
                    density[j, i, s2, s1] = sum(exp.(-beta * Es) .* wfs[:, s1] .* wfs[:, s2]) / Z
                end
            end
        end
    end

    DiagonalDensity(density)
end

"""
Diagonal density for a small PIGS system.
"""
struct PigsDiagonalDensity <: AbstractDiagonalDensity
    "Diagonal density."
    density::Array{Float64,4}
    "Exact diagonal density."
    density_exact::Array{Float64,4}
end

"""
    PigsDiagonalDensity(sys::System{S,2}, trial::TrialWavefunction{S,2}, beta::Float64, basis_size::Int, extent::NTuple{4,Float64}; lengths::NTuple{2,Int}=(101, 101), progress_output::IO=stderr)

Compute the diagonal (in position) density matrix for the 2-mode system `sys`
with `trial` propagated by `beta` using `basis_size` basis functions.

The extent for each mode is given in `extent`, and the number of points is
given in `lengths`.

The progress meter is written to `progress_output`.
"""
function PigsDiagonalDensity(sys::System{S,2}, trial::TrialWavefunction{S,2}, beta::Float64, basis_size::Int, extent::NTuple{4,Float64}; lengths::NTuple{2,Int}=(101, 101), progress_output::IO=stderr) where {S}
    basis = Basis(sys, basis_size)
    h0, V = operators(basis, sys)
    basis_vectors = vectors(basis)

    F = eigen(Symmetric(h0 + V))
    Es = F.values
    Vs = F.vectors

    # Trial wavefunction.
    trial_H = Vs' * trial_mode(trial, basis)

    # Propagated wavefunction.
    wf_H = exp.(-0.5 * beta * Es) .* trial_H
    wf_vec = Vs * wf_H

    # Exact wavefunction.
    wf_exact = Vs[:, 1]

    Z = dot(wf_H, wf_H)

    M = 2
    q1_min, q1_max, q2_min, q2_max = extent
    q1_length, q2_length = lengths
    density = zeros(q2_length, q1_length, S, S)
    density_exact = zeros(q2_length, q1_length, S, S)

    meter = Progress(q1_length, output=progress_output)
    for (i, q1) in ProgressWrapper(enumerate(range(q1_min; stop=q1_max, length=q1_length)), meter)
        for (j, q2) in enumerate(range(q2_min; stop=q2_max, length=q2_length))
            wfs = zeros(Float64, S)
            wfs_exact = zeros(Float64, S)
            for idx in 1:basis.dim
                s = basis_vectors[M+1, idx]
                ho_wf_eval = ho_wf(basis_vectors[1:M, idx], [q1, q2])
                wfs[s] += ho_wf_eval * wf_vec[idx]
                wfs_exact[s] += ho_wf_eval * wf_exact[idx]
            end
            density[j, i, :, :] = wfs * wfs' / Z
            density_exact[j, i, :, :] = wfs_exact * wfs_exact'
        end
    end

    PigsDiagonalDensity(density, density_exact)
end
