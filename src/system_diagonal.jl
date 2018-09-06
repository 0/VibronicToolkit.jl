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

"""
    DiagonalSystem(path::AbstractString)

Create a diagonal system from the JSON description in the file at `path`.

If there is any coupling between surfaces, `SurfaceCouplingException` is
thrown.
"""
DiagonalSystem(path::AbstractString) = path |> DenseSystem |> DiagonalSystem

diag(sys::DiagonalSystem{S,M}) where {S,M} = sys

isdiag(sys::DiagonalSystem{S,M}) where {S,M} = true

function JSON.lower(sys::DiagonalSystem{S,M}) where {S,M}
    result = Dict()

    result["number of surfaces"] = S
    result["number of modes"] = M

    result["energies"] = diag(sys.energy)
    result["frequencies"] = sys.freq[:, 1]
    result["linear couplings"] = permutedims(diag(sys.lin, 2, 3))
    result["quadratic couplings"] = permutedims(diag(sys.quad, 3, 4), [3, 2, 1])

    result
end
