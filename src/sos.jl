# Sum over states solution for systems with few modes.

"""
Sum over states solution for a small system.
"""
abstract type AbstractSumOverStates <: Solution end

"""
Sum over states solution for a small system at finite temperature.
"""
struct SumOverStates <: AbstractSumOverStates
    "Partition function."
    Z::Float64
    "Energy."
    E::Float64
    "Heat capacity."
    Cv::Float64
end

"""
    SumOverStates(sys::System, beta::Float64, basis_size::Int)

Calculate the solution for `sys` at `beta` using `basis_size` basis functions.
"""
function SumOverStates(sys::System, beta::Float64, basis_size::Int)
    basis = Basis(sys, basis_size)
    h0, V = operators(basis, sys)

    Es = eigvals(Symmetric(h0 + V))

    Z = sum(exp.(-beta * Es))
    E = sum(exp.(-beta * Es) .* Es) / Z
    Cv = sum(exp.(-beta * Es) .* (Es .- E).^2) / Z * beta^2

    SumOverStates(Z, E, Cv)
end

"""
Sum over states solution for a small PIGS system.
"""
struct PigsSumOverStates <: AbstractSumOverStates
    "Pseudo-partition function."
    Z::Float64
    "Energy."
    E::Float64
    "Von Neumann entanglement entropy."
    SvN::Float64
    "Order-2 Rényi entanglement entropy."
    S2::Float64

    "Exact ground state energy."
    E0_exact::Float64
    "Exact ground state von Neumann entanglement entropy."
    SvN_exact::Float64
    "Exact ground state order-2 Rényi entanglement entropy."
    S2_exact::Float64
end

"""
    PigsSumOverStates(sys::System{S,M}, trial::TrialWavefunction{S,M}, beta::Float64, basis_size::Int)

Calculate the solution for `sys` with `trial` propagated by `beta` using
`basis_size` basis functions.
"""
function PigsSumOverStates(sys::System{S,M}, trial::TrialWavefunction{S,M}, beta::Float64, basis_size::Int) where {S,M}
    basis = Basis(sys, basis_size)
    h0, V = operators(basis, sys)

    F = eigen(Symmetric(h0 + V))
    Es = F.values
    Vs = F.vectors

    # Trial wavefunction.
    trial_H = Vs' * trial_mode(trial, basis)

    # Propagated wavefunction.
    wf_H = exp.(-0.5 * beta * Es) .* trial_H
    wf_vec = Vs * wf_H

    # Density matrix and reduced density matrix.
    rho = wf_vec * wf_vec'
    rho_surface = ptrace_modes(basis, rho)

    Z = dot(wf_H, wf_H)
    E = dot(wf_H, Es .* wf_H) / Z
    SvN = S_vn(rho_surface / Z)
    S2 = S_renyi(rho_surface / Z)

    # Exact wavefunction, density matrix, and reduced density matrix.
    wf_exact = Vs[:, 1]
    rho_exact = wf_exact * wf_exact'
    rho_surface_exact = ptrace_modes(basis, rho_exact)

    E0_exact = Es[1]
    SvN_exact = S_vn(rho_surface_exact)
    S2_exact = S_renyi(rho_surface_exact)

    PigsSumOverStates(Z, E, SvN, S2, E0_exact, SvN_exact, S2_exact)
end
