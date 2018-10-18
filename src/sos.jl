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

Calculate the solution for `sys` at `beta` with `basis_size` basis functions.
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

    "Exact ground state energy."
    E0_exact::Float64
end

"""
    PigsSumOverStates(sys::System, beta::Float64, basis_size::Int)

Calculate the solution for `sys` at `beta` with `basis_size` basis functions
and a uniform (in space and surfaces) trial wavefunction.
"""
function PigsSumOverStates(sys::System, beta::Float64, basis_size::Int)
    basis = Basis(sys, basis_size)
    h0, V = operators(basis, sys)

    F = eigen(Symmetric(h0 + V))
    Es = F.values
    Vs = F.vectors

    # Uniform trial wavefunction represented in the eigenbasis of the
    # Hamiltonian.
    trial = Vs' * trial_uniform(basis)

    Z = sum(exp.(-beta * Es) .* abs2.(trial))
    E = sum(exp.(-beta * Es) .* Es .* abs2.(trial)) / Z
    E0_exact = Es[1]

    PigsSumOverStates(Z, E, E0_exact)
end
