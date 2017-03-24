# Operators in the harmonic oscillator basis.

"""
Harmonic oscillator basis for `S` surfaces and `M` modes.
"""
immutable Basis{S,M}
    "Single-mode basis size."
    size::Int

    "Single-surface basis size."
    dim1::Int
    "Total basis size."
    dim::Int
end

"""
    Basis{S,M}(::System{S,M}, size::Int)

Create a basis with `size` basis functions for a single mode.
"""
function Basis{S,M}(::System{S,M}, size::Int)
    # At least 1 basis function.
    size >= 1 || throw(DomainError())

    dim1 = size^M
    dim = S * dim1

    Basis{S,M}(size, dim1, dim)
end

"""
    id(basis::Basis)

Identity operator in `basis` for a single mode.
"""
id(basis::Basis) = eye(basis.size)

"""
    a(basis::Basis)

Annihilation operator in `basis` for a single mode.
"""
a(basis::Basis) = diagm(sqrt(1.0:basis.size-1), 1)

"""
    n(basis::Basis)

Number operator in `basis` for a single mode.
"""
n(basis::Basis) = diagm(0.0:basis.size-1)

"""
    q(basis::Basis)

Position operator in `basis` for a single mode.
"""
q(basis::Basis) = (a(basis)' + a(basis)) / sqrt(2.0)

"""
    mkop{S,M}(basis::Basis{S,M}, op::Array{Float64,2}, idx::Int)

Create a multi-mode operator in `basis` that is the tensor product of `op` at
mode `idx` and the identity operator at all other modes.
"""
function mkop{S,M}(basis::Basis{S,M}, op::Array{Float64,2}, idx::Int)
    # Nothing to do for a single mode.
    M == 1 && return op

    ops = []
    for m in 1:M
        if m == idx
            push!(ops, op)
        else
            push!(ops, id(basis))
        end
    end

    result = kron(ops...)
    size(result) == (basis.dim1, basis.dim1) || warn("Bad size in mkop")

    result
end

"""
    mkop{S,M}(basis::Basis{S,M}, op, idx::Int)

Create a multi-mode operator in `basis` that is the tensor product of `op` at
mode `idx` and the identity operator at all other modes, where `op` is a
function that generates the operator in the basis.
"""
mkop{S,M}(basis::Basis{S,M}, op, idx::Int) = mkop(basis, op(basis), idx)

"""
    mkid{S,M}(basis::Basis{S,M})

Create a multi-mode identity operator in `basis`.
"""
mkid{S,M}(basis::Basis{S,M}) = mkop(basis, id(basis), 1)

"""
    operators{S,M}(basis::Basis{S,M}, sys::System{S,M})

Create `h0` and `V` operators in `basis` for the system described by `sys`.
"""
function operators{S,M}(basis::Basis{S,M}, sys::System{S,M})
    # It's simpler to populate the values into higher rank tensors.
    h0s = zeros(basis.dim1, basis.dim1, S, S)
    Vs = zeros(basis.dim1, basis.dim1, S, S)

    for s1 in 1:S
        for s2 in 1:S
            if s1 == s2
                h0s[:, :, s1, s1] += sys.energy[s1, s1] * mkid(basis)

                for m in 1:M
                    h0s[:, :, s1, s1] += sys.freq[m, s1] * (mkop(basis, n, m) + 0.5 * mkop(basis, id, m))
                end
            else
                Vs[:, :, s2, s1] += sys.energy[s2, s1] * mkid(basis)
            end

            for m in 1:M
                if s1 == s2
                    h0s[:, :, s1, s1] += sys.lin[m, s1, s1] * mkop(basis, q, m)
                else
                    Vs[:, :, s2, s1] += sys.lin[m, s2, s1] * mkop(basis, q, m)
                end
            end

            for m1 in 1:M
                for m2 in 1:M
                    Vs[:, :, s2, s1] += 0.5 * sys.quad[m2, m1, s2, s1] * mkop(basis, q, m2) * mkop(basis, q, m1)
                end
            end
        end
    end

    # Flatten into matrices.
    h0 = reshape(permutedims(h0s, [1, 3, 2, 4]), (basis.dim, basis.dim))
    maximum(abs(h0' - h0)) < 1e-13 || warn("Asymmetric h0: $(maximum(abs(h0' - h0)))")
    V = reshape(permutedims(Vs, [1, 3, 2, 4]), (basis.dim, basis.dim))
    maximum(abs(V' - V)) < 1e-13 || warn("Asymmetric V: $(maximum(abs(V' - V)))")

    h0, V
end
