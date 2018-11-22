# Description of a vibronic system that cannot have inter-surface coupling.

"""
System of `M` coupled harmonic oscillators (modes) across `S` surfaces, with no
coupling between surfaces.
"""
struct DiagonalSystem{S,M} <: System{S,M}
    "Frequencies (M, S)."
    freq::Matrix{Float64}
    "Hamiltonian coefficients."
    coef::HamiltonianCoefficients{S,M}

    "Energy offsets due to linear terms (S)."
    deltas::Vector{Float64}
    "Position offsets due to linear terms (M, S)."
    ds::Matrix{Float64}

    function DiagonalSystem(freq::AbstractMatrix{Float64}, coef::HamiltonianCoefficients{S,M}) where {S,M}
        size(freq) == (M, S) || throw(DomainError(size(freq), "Unexpected freq dimensions."))
        isdiag(coef) || throw(SurfaceCouplingException())

        deltas = zeros(S)
        ds = zeros(M, S)
        for s in 1:S
            for m in 1:M
                deltas[s] += -0.5 * coef[1][m, s, s].^2 / freq[m, s]
                ds[m, s] = -coef[1][m, s, s] / freq[m, s]
            end
        end

        new{S,M}(freq, coef, deltas, ds)
    end
end

function DiagonalSystem(freq::AbstractMatrix{Float64}, coef::AbstractVector)
    (M, S) = size(freq)
    DiagonalSystem(freq, HamiltonianCoefficients{S,M}(coef...))
end

bare(::Type{T}) where {T<:DiagonalSystem} = DiagonalSystem

diag(sys::DiagonalSystem) = sys

isdiag(sys::DiagonalSystem) = true
