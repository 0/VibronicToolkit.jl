# Trotter factorized solution for systems with few modes.

"""
Trotter factorized solution for a small system.
"""
immutable Trotter <: Solution
    "Partition function."
    Z::Float64
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
    A = expm(-tau * h0) * expm(-tau * V)

    Z = trace(A^P)

    Trotter(Z)
end
