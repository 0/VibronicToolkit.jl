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
    isdiag(sys::System{S,M})

Determine whether `sys` is diagonal (in content, rather than type or
structure). Returns `true` if there is no coupling between surfaces and `false`
otherwise.
"""
function isdiag end

include("system_dense.jl")
include("system_diagonal.jl")
include("system_io.jl")

"""
    DenseSystem(sys::DiagonalSystem{S,M})

Create a dense system from the diagonal system `sys`.
"""
function DenseSystem(sys::DiagonalSystem{S,M}) where {S,M}
    DenseSystem(sys.energy, sys.freq, sys.lin, sys.quad)
end

"""
    DiagonalSystem(sys::DenseSystem{S,M})

Create a diagonal system from the dense system `sys`.

If there is any coupling between surfaces, `SurfaceCouplingException` is
thrown.
"""
function DiagonalSystem(sys::DenseSystem{S,M}) where {S,M}
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
    simplify(sys::T)

Generate a simplified version of `sys` with no quadratic coupling.
"""
function simplify(sys::T) where {T<:System}
    T(sys.energy, sys.freq, sys.lin, zero(sys.quad))
end
