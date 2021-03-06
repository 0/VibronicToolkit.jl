# Sum over states solution for systems with few modes.

"""
Sum over states solution for a small system.
"""
struct SumOverStates <: Solution
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
