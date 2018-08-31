# Description of a vibronic system that may have inter-surface coupling.

"""
System of `M` coupled harmonic oscillators (modes) across `S` surfaces, with
the possibility for coupling between surfaces.
"""
struct DenseSystem{S,M} <: System{S,M}
    "Energy offsets (S, S)."
    energy::Matrix{Float64}
    "Frequencies (M, S)."
    freq::Matrix{Float64}
    "Linear prefactors (M, S, S)."
    lin::Array{Float64,3}
    "Quadratic prefactors (M, M, S, S)."
    quad::Array{Float64,4}

    function DenseSystem{S,M}(energy::AbstractMatrix{Float64}, freq::AbstractMatrix{Float64}, lin::AbstractArray{Float64,3}, quad::AbstractArray{Float64,4}) where {S,M}
        check_shape(S, M, energy, freq, lin, quad)
        new{S,M}(energy, freq, lin, quad)
    end

    function DenseSystem(energy::AbstractMatrix{Float64}, freq::AbstractMatrix{Float64}, lin::AbstractArray{Float64,3}, quad::AbstractArray{Float64,4})
        M, S = size(freq)
        DenseSystem{S,M}(energy, freq, lin, quad)
    end
end

function diag(sys::DenseSystem{S,M}) where {S,M}
    energy = Diagonal(sys.energy)
    freq = sys.freq
    lin = zero(sys.lin)
    for m in 1:M
        lin[m, :, :] .= Diagonal(sys.lin[m, :, :])
    end
    quad = zero(sys.quad)
    for m1 in 1:M
        for m2 in 1:M
            quad[m2, m1, :, :] .= Diagonal(sys.quad[m2, m1, :, :])
        end
    end
    DiagonalSystem(energy, freq, lin, quad)
end

isdiag(sys::DenseSystem{S,M}) where {S,M} = isdiag(sys.energy, sys.lin, sys.quad)
