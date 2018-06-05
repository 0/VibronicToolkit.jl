# Description of the vibronic system.

"""
System of `M` coupled harmonic oscillators (modes) across `S` surfaces.
"""
struct System{S,M}
    "Energy offsets (S, S)."
    energy::Matrix{Float64}
    "Frequencies (M, S)."
    freq::Matrix{Float64}
    "Linear prefactors (M, S, S)."
    lin::Array{Float64,3}
    "Quadratic prefactors (M, M, S, S)."
    quad::Array{Float64,4}
end

"""
    System(path)

Create a system from the JSON description in the file at `path`.
"""
function System(path)
    data = JSON.parsefile(path)

    S = data["number of surfaces"]
    # At least 1 surface.
    S >= 1 || throw(DomainError())

    M = data["number of modes"]
    # At least 1 mode.
    M >= 1 || throw(DomainError())

    energy = zeros(S, S)
    if haskey(data, "energies")
        length(data["energies"]) == size(energy, 1) || error("Bad energy size")
        for (idx1, data1) in enumerate(data["energies"])
            length(data1) == size(energy, 2) || error("Bad energy size")
            for (idx2, data2) in enumerate(data1)
                energy[idx1, idx2] = data2
            end
        end
    end
    all(isapprox.(energy, energy'; rtol=1e-12)) || error("Asymmetric energy")

    freq = zeros(M, S)
    # We repeat the frequency values over the surfaces.
    if haskey(data, "frequencies")
        length(data["frequencies"]) == size(freq, 1) || error("Bad freq size")
        for (idx1, data1) in enumerate(data["frequencies"])
            freq[idx1, :] = data1
        end
    end
    all(freq .> 0) || error("Nonpositive freq")

    lin = zeros(M, S, S)
    if haskey(data, "linear couplings")
        length(data["linear couplings"]) == size(lin, 1) || error("Bad lin size")
        for (idx1, data1) in enumerate(data["linear couplings"])
            length(data1) == size(lin, 2) || error("Bad lin size")
            for (idx2, data2) in enumerate(data1)
                length(data2) == size(lin, 3) || error("Bad lin size")
                for (idx3, data3) in enumerate(data2)
                    lin[idx1, idx2, idx3] = data3
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
                    length(data3) == size(quad, 4) || error("Bad quad size")
                    for (idx4, data4) in enumerate(data3)
                        quad[idx1, idx2, idx3, idx4] = data4
                    end
                end
            end
        end
    end
    all(isapprox.(quad, permutedims(quad, [1, 2, 4, 3]); rtol=1e-12)) || error("Asymmetric quad")

    System{S,M}(energy, freq, lin, quad)
end

function JSON.lower{S,M}(sys::System{S,M})
    result = Dict()

    result["number of surfaces"] = S
    result["number of modes"] = M

    result["energies"] = sys.energy'
    result["frequencies"] = sys.freq[:, 1]
    result["linear couplings"] = permutedims(sys.lin, [3, 2, 1])
    result["quadratic couplings"] = permutedims(sys.quad, [4, 3, 2, 1])

    result
end

Base.write{S,M}(io::IO, sys::System{S,M}) = JSON.print(io, sys)

"""
    simplify{S,M}(sys::System{S,M})

Generate a simplified version of `sys` with no quadratic coupling and no
inter-surface linear coupling.
"""
function simplify{S,M}(sys::System{S,M})
    energy_new = diagm(diag(sys.energy))

    lin_new = zeros(sys.lin)
    for s in 1:S
        lin_new[:, s, s] .= sys.lin[:, s, s]
    end

    System{S,M}(energy_new, sys.freq, lin_new, zeros(sys.quad))
end

"""
    is_coupled{S,M}(sys::System{S,M})

Whether `sys` has coupling between surfaces.
"""
function is_coupled{S,M}(sys::System{S,M})
    isdiag(sys.energy) || return true

    for m in 1:M
        isdiag(sys.lin[m, :, :]) || return true
    end

    for m1 in 1:M
        for m2 in 1:M
            isdiag(sys.quad[m2, m1, :, :]) || return true
        end
    end

    false
end
