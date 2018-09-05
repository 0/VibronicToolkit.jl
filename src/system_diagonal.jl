# Description of a vibronic system that cannot have inter-surface coupling.

"""
System of `M` coupled harmonic oscillators (modes) across `S` surfaces, with no
coupling between surfaces.
"""
struct DiagonalSystem{S,M} <: System{S,M}
    "Energy offsets (S, S)."
    energy::Matrix{Float64}
    "Frequencies (M, S)."
    freq::Matrix{Float64}
    "Linear prefactors (M, S, S)."
    lin::Array{Float64,3}
    "Quadratic prefactors (M, M, S, S)."
    quad::Array{Float64,4}

    "Energy offsets due to linear terms (S)."
    deltas::Vector{Float64}
    "Position offsets due to linear terms (M, S)."
    ds::Matrix{Float64}

    function DiagonalSystem{S,M}(energy::AbstractMatrix{Float64}, freq::AbstractMatrix{Float64}, lin::AbstractArray{Float64,3}, quad::AbstractArray{Float64,4}) where {S,M}
        check_shape(S, M, energy, freq, lin, quad)
        isdiag(energy, lin, quad) || throw(SurfaceCouplingException())

        deltas = zeros(S)
        ds = zeros(M, S)
        for s in 1:S
            for m in 1:M
                deltas[s] += -0.5 * lin[m, s, s].^2 / freq[m, s]
                ds[m, s] = -lin[m, s, s] / freq[m, s]
            end
        end

        new{S,M}(energy, freq, lin, quad, deltas, ds)
    end

    function DiagonalSystem(energy::AbstractMatrix{Float64}, freq::AbstractMatrix{Float64}, lin::AbstractArray{Float64,3}, quad::AbstractArray{Float64,4})
        M, S = size(freq)
        DiagonalSystem{S,M}(energy, freq, lin, quad)
    end
end

diag(sys::DiagonalSystem) = sys

isdiag(sys::DiagonalSystem) = true
