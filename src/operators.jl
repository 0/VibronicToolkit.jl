# Operators in the harmonic oscillator basis.

"""
Harmonic oscillator basis for `S` surfaces and `M` modes.
"""
struct Basis{S,M}
    "Single-mode basis size."
    size::Int

    "Single-surface basis size."
    dim1::Int
    "Total basis size."
    dim::Int
end

"""
    Basis(::System, size::Int)

Create a basis with `size` basis functions for a single mode.
"""
function Basis(::System{S,M}, size::Int) where {S,M}
    size >= 1 || throw(DomainError(size, "At least 1 basis function."))

    dim1 = size^M
    dim = S * dim1

    Basis{S,M}(size, dim1, dim)
end

"""
    id(basis::Basis)

Identity operator in `basis` for a single mode.
"""
id(basis::Basis) = Matrix{Float64}(I, basis.size, basis.size)

"""
    a(basis::Basis)

Annihilation operator in `basis` for a single mode.
"""
a(basis::Basis) = diagm(1 => sqrt.(1.0:basis.size-1))

"""
    n(basis::Basis)

Number operator in `basis` for a single mode.
"""
n(basis::Basis) = diagm(0 => 0.0:basis.size-1)

"""
    q(basis::Basis)

Position operator in `basis` for a single mode.
"""
q(basis::Basis) = (a(basis)' + a(basis)) / sqrt(2.0)

"""
    mkop(basis::Basis, op::Matrix{Float64}, idx::Int)

Create a multi-mode operator in `basis` that is the tensor product of `op` at
mode `idx` and the identity operator at all other modes.
"""
function mkop(basis::Basis{S,M}, op::Matrix{Float64}, idx::Int) where {S,M}
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
    size(result) == (basis.dim1, basis.dim1) || @warn "Bad size in mkop"

    result
end

"""
    mkop(basis::Basis, op, idx::Int)

Create a multi-mode operator in `basis` that is the tensor product of `op` at
mode `idx` and the identity operator at all other modes, where `op` is a
function that generates the operator in the basis.
"""
mkop(basis::Basis, op, idx::Int) = mkop(basis, op(basis), idx)

"""
    mkid(basis::Basis)

Create a multi-mode identity operator in `basis`.
"""
mkid(basis::Basis) = mkop(basis, id(basis), 1)

"""
    operators(basis::Basis{S,M}, sys::System{S,M}; splitting::Splitting=h0_V)

Split the Hamiltonian `H` for the system described by `sys` into operators `A`
and `B` in `basis` so that `H = A + B`. The choice of operators is determined
by `splitting`.
"""
function operators(basis::Basis{S,M}, sys::System{S,M}; splitting::Splitting=h0_V) where {S,M}
    # It's simpler to populate the values into higher rank tensors.
    As = zeros(basis.dim1, basis.dim1, S, S)
    Bs = zeros(basis.dim1, basis.dim1, S, S)

    for s in 1:S
        for m in 1:M
            As[:, :, s, s] .+= sys.freq[m, s] * (mkop(basis, n, m) + 0.5 * mkop(basis, id, m))

            if splitting == p2_U
                As[:, :, s, s] .-= sys.freq[m, s] * 0.5 * mkop(basis, q, m)^2
                Bs[:, :, s, s] .+= sys.freq[m, s] * 0.5 * mkop(basis, q, m)^2
            end
        end
    end

    coefficients(sys.coef) do ord, s1, s2, idx, val
        k = val * mkid(basis)
        for m in idx
            k *= mkop(basis, q, m)
        end
        if ord <= 1 && s1 == s2
            if splitting == h0_V
                As[:, :, s1, s1] .+= k
            elseif splitting == p2_U
                Bs[:, :, s1, s1] .+= k
            else
                error("Unsupported splitting: $(splitting)")
            end
        else
            Bs[:, :, s2, s1] .+= k
        end
    end

    # Flatten into matrices.
    A = reshape(permutedims(As, [1, 3, 2, 4]), (basis.dim, basis.dim))
    maximum(abs.(A' - A)) < 1e-13 || @warn "Asymmetric A: $(maximum(abs.(A' - A)))"
    B = reshape(permutedims(Bs, [1, 3, 2, 4]), (basis.dim, basis.dim))
    maximum(abs.(B' - B)) < 1e-13 || @warn "Asymmetric B: $(maximum(abs.(B' - B)))"

    A, B
end

"""
    vectors(basis::Basis{S,M})

Generate an array of basis vectors for `basis` in the order corresponding to
the layout of the operator matrices.

Each column is a basis vector, with the first `M` values giving the (0-indexed)
state labels of the modes, and the final value giving the (1-indexed) surface
label.
"""
function vectors(basis::Basis{S,M}) where {S,M}
    result = Array{Int}(undef, M+1, basis.dim)

    conf = zeros(Int, M)
    idx = 1
    for s in 1:S
        first = true
        while true
            for m in M:-1:1
                !first && (conf[m] += 1)
                if conf[m] < basis.size
                    break
                else
                    conf[m] = 0
                end
            end
            if first
                first = false
            elseif iszero(conf)
                break
            end

            result[1:M, idx] .= conf
            result[M+1, idx] = s

            idx += 1
        end
    end
    idx-1 == basis.dim || @warn "Invalid number of basis vectors: $(idx-1)"

    result
end

"""
    ptrace_modes(basis::Basis, rho::AbstractMatrix{Float64})

Perform the partial trace over the modes for `rho` in `basis`.
"""
ptrace_modes(basis::Basis{S,M}, rho::AbstractMatrix{Float64}) where {S,M} = ptrace(rho, (basis.dim1, S), 1)
