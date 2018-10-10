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

include("system_dense.jl")
include("system_diagonal.jl")
include("system_io.jl")

function Base.:(==)(sys1::System{S,M}, sys2::System{S,M}) where {S,M}
    sys1.energy == sys2.energy || return false
    sys1.freq == sys2.freq || return false
    sys1.lin == sys2.lin || return false
    sys1.quad == sys2.quad || return false
    true
end

function Base.isapprox(sys1::System{S,M}, sys2::System{S,M}; kwargs...) where {S,M}
    isapprox(sys1.energy, sys2.energy; kwargs...) || return false
    isapprox(sys1.freq, sys2.freq; kwargs...) || return false
    isapprox(sys1.lin, sys2.lin; kwargs...) || return false
    isapprox(sys1.quad, sys2.quad; kwargs...) || return false
    true
end

"""
    DenseSystem(sys::DiagonalSystem)

Create a dense system from the diagonal system `sys`.
"""
DenseSystem(sys::DiagonalSystem) = DenseSystem(sys.energy, sys.freq, sys.lin, sys.quad)

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

"""
    isdiag(energy::AbstractMatrix{Float64}, lin::AbstractArray{Float64,3}, quad::AbstractArray{Float64,4})

Determine whether the component tensors of a system are diagonal in surfaces.
"""
function isdiag(energy::AbstractMatrix{Float64}, lin::AbstractArray{Float64,3}, quad::AbstractArray{Float64,4})
    isdiag(energy) || return false
    for m in 1:size(lin, 1)
        isdiag(lin[m, :, :]) || return false
    end
    for m1 in 1:size(quad, 2)
        for m2 in 1:size(quad, 1)
            isdiag(quad[m2, m1, :, :]) || return false
        end
    end
    true
end

"""
    check_shape(S::Int, M::Int, energy::AbstractMatrix{Float64}, freq::AbstractMatrix{Float64}, lin::AbstractArray{Float64,3}, quad::AbstractArray{Float64,4})

Throw the appropriate exception if the component tensors of a system have the
wrong shape.
"""
function check_shape(S::Int, M::Int, energy::AbstractMatrix{Float64}, freq::AbstractMatrix{Float64}, lin::AbstractArray{Float64,3}, quad::AbstractArray{Float64,4})
    S >= 1 || throw(DomainError(S, "At least 1 surface."))
    M >= 1 || throw(DomainError(M, "At least 1 mode."))
    size(energy) == (S, S) || throw(DomainError(size(energy), "Expected energy dimensions: $(S), $(S)."))
    size(freq) == (M, S) || throw(DomainError(size(freq), "Expected freq dimensions: $(M), $(S)."))
    size(lin) == (M, S, S) || throw(DomainError(size(lin), "Expected lin dimensions: $(M), $(S), $(S)."))
    size(quad) == (M, M, S, S) || throw(DomainError(size(quad), "Expected quad dimensions: $(M), $(M), $(S), $(S)."))
    nothing
end

"""
    simplify(sys::System)

Generate a simplified version of `sys` with no quadratic coupling.
"""
function simplify(sys::T) where {T<:System}
    T(sys.energy, sys.freq, sys.lin, zero(sys.quad))
end

"""
    issimple(sys::System)

Determine whether `sys` is simple (has no quadratic coupling).
"""
issimple(sys::System) = iszero(sys.quad)

"""
    potential(sys::System, qs::Vector{Float64})

Compute the potential matrix (full Hamiltonian without kinetic energy) for
`sys` at `qs`.
"""
function potential(sys::System{S,M}, qs::Vector{Float64}) where {S,M}
    V = copy(sys.energy)
    for s1 in 1:S
        for s2 in 1:S
            for m in 1:M
                if s1 == s2
                    V[s2, s1] += 0.5 * sys.freq[m] * qs[m]^2
                end
                V[s2, s1] += sys.lin[m, s2, s1] * qs[m]
            end
            for m1 in 1:M
                for m2 in 1:M
                    V[s2, s1] += 0.5 * sys.quad[m2, m1, s2, s1] * qs[m2] * qs[m1]
                end
            end
        end
    end
    Symmetric(V)
end
