# Trotter factorized solution for systems with few modes.

"""
Trotter factorized solution for a small system.
"""
immutable Trotter <: Solution
    "Partition function."
    Z::Float64
    "Energy."
    E::Float64
    "Heat capacity."
    Cv::Float64
end

"""
    Trotter{S,M}(sys::System{S,M}, beta::Beta, basis_size::Int, P::Int)

Calculate the solution for `sys` at `beta` with `basis_size` basis functions
and `P` links.
"""
function Trotter{S,M}(sys::System{S,M}, beta::Beta, basis_size::Int, P::Int)
    basis = Basis(sys, basis_size)
    h0, V = operators(basis, sys)

    # Single Trotter product.
    tau = beta / P
    rho = expm(-tau * h0) * expm(-tau * V)

    # Full path.
    path = rho^P

    Z = trace(path)
    E = trace(path * (h0 + V)) / Z
    Cv = (trace(path * (h0 + V)^2) / Z - E^2) * beta^2

    Trotter(Z, E, Cv)
end
