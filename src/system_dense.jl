# Description of a vibronic system that may have inter-surface coupling.

"""
System of `M` coupled harmonic oscillators (modes) across `S` surfaces, with
the possibility for coupling between surfaces.
"""
struct DenseSystem{S,M} <: System{S,M}
    "Frequencies (M, S)."
    freq::Matrix{Float64}
    "Hamiltonian coefficients."
    coef::HamiltonianCoefficients{S,M}

    function DenseSystem(freq::AbstractMatrix{Float64}, coef::HamiltonianCoefficients{S,M}) where {S,M}
        size(freq) == (M, S) || throw(DomainError(size(freq), "Unexpected freq dimensions."))
        new{S,M}(freq, coef)
    end
end

function DenseSystem(freq::AbstractMatrix{Float64}, coef::AbstractVector)
    (M, S) = size(freq)
    DenseSystem(freq, HamiltonianCoefficients{S,M}(coef...))
end

function DenseSystem(freq::AbstractMatrix{Float64}, coef)
    (M, S) = size(freq)
    DenseSystem(freq, HamiltonianCoefficients{S,M}(coef))
end

bare(::Type{T}) where {T<:DenseSystem} = DenseSystem

function diag(sys::DenseSystem{S,M}) where {S,M}
    coef_d = Dict{Int,Array{Float64}}()
    for (ord, val) in sys.coef
        val_diag = zero(val)
        for idx in mode_indices(val)
            val_diag[idx, :, :] .= Diagonal(val[idx, :, :])
        end
        coef_d[ord] = val_diag
    end
    DiagonalSystem(sys.freq, HamiltonianCoefficients{S,M}(coef_d))
end

isdiag(sys::DenseSystem) = isdiag(sys.coef)
