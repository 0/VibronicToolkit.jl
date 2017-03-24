# Sum over states solution for systems with few modes.

"""
Sum over states solution for a small system.
"""
immutable SumOverStates <: Solution
    "Partition function."
    Z::Float64
    "Energy."
    E::Float64
    "Heat capacity."
    Cv::Float64
end

"""
    SumOverStates{S,M}(sys::System{S,M}, beta::Beta, basis_size::Int)

Calculate the solution for `sys` at `beta` with `basis_size` basis functions.
"""
function SumOverStates{S,M}(sys::System{S,M}, beta::Beta, basis_size::Int)
    basis = Basis(sys, basis_size)
    h0, V = operators(basis, sys)

    eigen = eigfact(Symmetric(h0 + V))
    Es = eigen[:values]

    Z = sum(exp(-beta * Es))
    E = sum(exp(-beta * Es) .* Es) / Z
    Cv = sum(exp(-beta * Es) .* (Es - E).^2) / Z

    SumOverStates(Z, E, Cv)
end
