# Description of a vibronic system.

"""
The system has non-zero off-diagonal elements that couple the surfaces.
"""
struct SurfaceCouplingException <: Exception end

"""
System of `M` coupled harmonic oscillators (modes) across `S` surfaces.
"""
abstract type System{S,M} end

"""
    diag(sys::System{S,M})

Extract the diagonal part of `sys` as a `DiagonalSystem{S,M}`.
"""
function diag end

"""
    isdiag(sys::System)

Determine whether `sys` is diagonal (in content, rather than type or
structure). Returns `true` if there is no coupling between surfaces and `false`
otherwise.
"""
function isdiag end

include("system_coef.jl")
include("system_dense.jl")
include("system_diagonal.jl")
include("system_io.jl")

function Base.:(==)(sys1::System{S,M}, sys2::System{S,M}) where {S,M}
    sys1.freq == sys2.freq || return false
    keys(sys1.coef) == keys(sys2.coef) || return false
    for ((ord1, val1), (ord2, val2)) in zip(sys1.coef, sys2.coef)
        ord1 == ord2 || throw(ErrorException("Order mismatch."))
        val1 == val2 || return false
    end
    true
end

function Base.isapprox(sys1::System{S,M}, sys2::System{S,M}; kwargs...) where {S,M}
    isapprox(sys1.freq, sys2.freq; kwargs...) || return false
    keys(sys1.coef) == keys(sys2.coef) || return false
    for ((ord1, val1), (ord2, val2)) in zip(sys1.coef, sys2.coef)
        ord1 == ord2 || throw(ErrorException("Order mismatch."))
        isapprox(val1, val2; kwargs...) || return false
    end
    true
end

function Base.getproperty(sys::System, sym::Symbol)
    if sym === :energy
        return sys.coef[0]
    elseif sym === :lin
        return sys.coef[1]
    elseif sym === :quad
        return sys.coef[2]
    else
        return getfield(sys, sym)
    end
end

"""
    DenseSystem(sys::DiagonalSystem)

Create a dense system from the diagonal system `sys`.
"""
DenseSystem(sys::DiagonalSystem) = DenseSystem(sys.freq, sys.coef)

"""
    DiagonalSystem(sys::DenseSystem)

Create a diagonal system from the dense system `sys`.

If there is any coupling between surfaces, `SurfaceCouplingException` is
thrown.
"""
function DiagonalSystem(sys::DenseSystem)
    isdiag(sys) || throw(SurfaceCouplingException())
    diag(sys)
end

Base.:*(x::Float64, sys::T) where {T<:System} = bare(T)(sys.freq, x * sys.coef)

"""
    simplify(sys::System; ord::Int=1)

Generate a simplified version of `sys` with no coefficients beyond `ord`.
"""
function simplify(sys::T; ord::Int=1) where {T<:System{S,M}} where {S,M}
    d = Dict{Int,Array{Float64}}(x for x in sys.coef if x.first <= ord)
    bare(T)(sys.freq, HamiltonianCoefficients{S,M}(d))
end

"""
    issimple(sys::System; ord::Int=1)

Determine whether `sys` is simple (has no coefficients beyond `ord`).
"""
function issimple(sys::System; ord::Int=1)
    ks = keys(sys.coef)
    isempty(ks) && return true
    maximum(ks) <= ord
end

"""
    potential(sys::System, qs::Vector{Float64})

Compute the potential matrix (full Hamiltonian without kinetic energy) for
`sys` at `qs`.
"""
function potential(sys::System{S,M}, qs::Vector{Float64}) where {S,M}
    V = zeros(Float64, S, S)
    for s in 1:S
            for m in 1:M
                    V[s, s] += 0.5 * sys.freq[m, s] * qs[m]^2
            end
    end
    coefficients(sys.coef) do ord, s1, s2, idx, val
        for m in idx
            val *= qs[m]
        end
        V[s2, s1] += val
    end
    Symmetric(V)
end
