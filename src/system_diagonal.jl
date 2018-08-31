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

    function DiagonalSystem{S,M}(energy::AbstractMatrix{Float64}, freq::AbstractMatrix{Float64}, lin::AbstractArray{Float64,3}, quad::AbstractArray{Float64,4}) where {S,M}
        check_shape(S, M, energy, freq, lin, quad)
        isdiag(energy, lin, quad) || throw(SurfaceCouplingException())
        new{S,M}(energy, freq, lin, quad)
    end

    function DiagonalSystem(energy::AbstractMatrix{Float64}, freq::AbstractMatrix{Float64}, lin::AbstractArray{Float64,3}, quad::AbstractArray{Float64,4})
        M, S = size(freq)
        DiagonalSystem{S,M}(energy, freq, lin, quad)
    end
end

diag(sys::DiagonalSystem{S,M}) where {S,M} = sys

isdiag(sys::DiagonalSystem{S,M}) where {S,M} = true
