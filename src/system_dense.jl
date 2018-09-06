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

"""
    DenseSystem(path::AbstractString)

Create a dense system from the JSON description in the file at `path`.
"""
function DenseSystem(path::AbstractString)
    data = JSON.parsefile(path)

    S = data["number of surfaces"]
    S >= 1 || throw(DomainError(S, "At least 1 surface."))

    M = data["number of modes"]
    M >= 1 || throw(DomainError(M, "At least 1 mode."))

    for key in ["linear coupling", "quadratic coupling"]
        haskey(data, key) && error("System file contains outdated key '$(key)'.")
    end

    energy = zeros(S, S)
    if haskey(data, "energies")
        length(data["energies"]) == size(energy, 1) || error("Bad energy size")
        for (idx1, data1) in enumerate(data["energies"])
            if length(data1) == 1
                energy[idx1, idx1] = data1[1]
            elseif length(data1) == size(energy, 2)
                for (idx2, data2) in enumerate(data1)
                    energy[idx1, idx2] = data2
                end
            else
                error("Bad energy size")
            end
        end
    end
    all(isapprox.(energy, permutedims(energy); rtol=1e-12)) || error("Asymmetric energy")

    freq = zeros(M, S)
    # We repeat the frequency values over the surfaces.
    if haskey(data, "frequencies")
        length(data["frequencies"]) == size(freq, 1) || error("Bad freq size")
        for (idx1, data1) in enumerate(data["frequencies"])
            freq[idx1, :] .= data1
        end
    end
    all(freq .> 0) || error("Nonpositive freq")

    lin = zeros(M, S, S)
    if haskey(data, "linear couplings")
        length(data["linear couplings"]) == size(lin, 1) || error("Bad lin size")
        for (idx1, data1) in enumerate(data["linear couplings"])
            length(data1) == size(lin, 2) || error("Bad lin size")
            for (idx2, data2) in enumerate(data1)
                if length(data2) == 1
                    lin[idx1, idx2, idx2] = data2[1]
                elseif length(data2) == size(lin, 3)
                    for (idx3, data3) in enumerate(data2)
                        lin[idx1, idx2, idx3] = data3
                    end
                else
                    error("Bad lin size")
                end
            end
        end
    end
    all(isapprox.(lin, permutedims(lin, [1, 3, 2]); rtol=1e-12)) || error("Asymmetric lin")

    quad = zeros(M, M, S, S)
    if haskey(data, "quadratic couplings")
        length(data["quadratic couplings"]) == size(quad, 1) || error("Bad quad size")
        for (idx1, data1) in enumerate(data["quadratic couplings"])
            length(data1) == size(quad, 2) || error("Bad quad size")
            for (idx2, data2) in enumerate(data1)
                length(data2) == size(quad, 3) || error("Bad quad size")
                for (idx3, data3) in enumerate(data2)
                    if length(data3) == 1
                        quad[idx1, idx2, idx3, idx3] = data3[1]
                    elseif length(data3) == size(quad, 4)
                        for (idx4, data4) in enumerate(data3)
                            quad[idx1, idx2, idx3, idx4] = data4
                        end
                    else
                        error("Bad quad size")
                    end
                end
            end
        end
    end
    all(isapprox.(quad, permutedims(quad, [1, 2, 4, 3]); rtol=1e-12)) || error("Asymmetric quad")

    DenseSystem{S,M}(energy, freq, lin, quad)
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

function JSON.lower(sys::DenseSystem{S,M}) where {S,M}
    result = Dict()

    result["number of surfaces"] = S
    result["number of modes"] = M

    result["energies"] = permutedims(sys.energy)
    result["frequencies"] = sys.freq[:, 1]
    result["linear couplings"] = permutedims(sys.lin, [3, 2, 1])
    result["quadratic couplings"] = permutedims(sys.quad, [4, 3, 2, 1])

    result
end
