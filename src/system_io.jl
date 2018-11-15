# IO utilities for vibronic systems.

function Base.read(io::IO, ::Type{System})
    data = JSON.parse(io)

    known_keys = ["number of surfaces", "number of modes", "energies",
                  "frequencies", "linear couplings", "quadratic couplings"]
    unrecognized_keys = setdiff(keys(data), known_keys)
    isempty(unrecognized_keys) || @warn "Unrecognized keys: $(unrecognized_keys)."

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

    energy, freq, lin, quad
end

Base.read(io::IO, ::Type{DenseSystem}) = DenseSystem(read(io, System)...)

Base.read(io::IO, ::Type{DiagonalSystem}) = DiagonalSystem(read(io, System)...)

Base.write(io::IO, sys::System) = JSON.print(io, sys)

function JSON.lower(sys::DenseSystem{S,M}) where {S,M}
    all(sys.freq .== sys.freq[:, 1]) || throw(DomainError(sys.freq, "All surfaces must have the same frequencies."))

    result = Dict()

    result["number of surfaces"] = S
    result["number of modes"] = M

    result["energies"] = permutedims(sys.energy)
    result["frequencies"] = sys.freq[:, 1]
    result["linear couplings"] = permutedims(sys.lin, [3, 2, 1])
    result["quadratic couplings"] = permutedims(sys.quad, [4, 3, 2, 1])

    result
end

function JSON.lower(sys::DiagonalSystem{S,M}) where {S,M}
    all(sys.freq .== sys.freq[:, 1]) || throw(DomainError(sys.freq, "All surfaces must have the same frequencies."))

    result = Dict()

    result["number of surfaces"] = S
    result["number of modes"] = M

    result["energies"] = diag(sys.energy)
    result["frequencies"] = sys.freq[:, 1]
    result["linear couplings"] = permutedims(diag(sys.lin, 2, 3))
    result["quadratic couplings"] = permutedims(diag(sys.quad, 3, 4), [3, 2, 1])

    result
end

function show_tensors(io::IO, sys::System{S,M}) where {S,M}
    if !iszero(sys.energy)
        println(io, "energies:")
        show_matrix(io, sys.energy)
    end
    println(io, "frequencies:")
    for m in 1:M
        println(io, " mode $(m)")
        show_vector(io, sys.freq[m, :])
    end
    if !iszero(sys.lin)
        println(io, "linear couplings:")
        for m in 1:M
            iszero(sys.lin[m, :, :]) && continue
            println(io, " mode $(m)")
            show_matrix(io, sys.lin[m, :, :])
        end
    end
    if !iszero(sys.quad)
        println(io, "quadratic couplings:")
        for m1 in 1:M
            for m2 in 1:M
                iszero(sys.quad[m2, m1, :, :]) && continue
                println(io, " modes $(m1), $(m2)")
                show_matrix(io, sys.quad[m2, m1, :, :])
            end
        end
    end
    nothing
end

function Base.show(io::IO, sys::DenseSystem{S,M}) where {S,M}
    println(io, typeof(sys))
    println(io, "Dense system with $(S) surface$(S == 1 ? "" : "s") and $(M) mode$(M == 1 ? "" : "s").")
    if isdiag(sys)
        println(io, "Does not contain inter-surface coupling.")
    else
        println(io, "Contains inter-surface coupling.")
    end
    show_tensors(io, sys)
    nothing
end

function Base.show(io::IO, sys::DiagonalSystem{S,M}) where {S,M}
    println(io, typeof(sys))
    println(io, "Diagonal system with $(S) surface$(S == 1 ? "" : "s") and $(M) mode$(M == 1 ? "" : "s").")
    show_tensors(io, sys)
    println(io, "energy offsets:")
    show_vector(io, sys.deltas)
    println(io, "position offsets:")
    for m in 1:M
        println(io, " mode $(m)")
        show_vector(io, sys.ds[m, :])
    end
    nothing
end
