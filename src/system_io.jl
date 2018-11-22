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

    freq = zeros(M, S)
    # We repeat the frequency values over the surfaces.
    if haskey(data, "frequencies")
        length(data["frequencies"]) == size(freq, 1) || error("Bad freq size")
        for (idx1, data1) in enumerate(data["frequencies"])
            freq[idx1, :] .= data1
        end
    end
    all(freq .> 0) || error("Nonpositive freq")

    coef_d = Dict{Int,Array{Float64}}()

    if haskey(data, "energies")
        energy = zeros(S, S)
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
        all(isapprox.(energy, permutedims(energy); rtol=1e-12)) || error("Asymmetric energy")
        coef_d[0] = energy
    end

    if haskey(data, "linear couplings")
        lin = zeros(M, S, S)
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
        all(isapprox.(lin, permutedims(lin, [1, 3, 2]); rtol=1e-12)) || error("Asymmetric lin")
        coef_d[1] = lin
    end

    if haskey(data, "quadratic couplings")
        quad = zeros(M, M, S, S)
        length(data["quadratic couplings"]) == size(quad, 1) || error("Bad quad size")
        for (idx1, data1) in enumerate(data["quadratic couplings"])
            length(data1) == size(quad, 2) || error("Bad quad size")
            for (idx2, data2) in enumerate(data1)
                length(data2) == size(quad, 3) || error("Bad quad size")
                for (idx3, data3) in enumerate(data2)
                    if length(data3) == 1
                        quad[idx1, idx2, idx3, idx3] = 0.5 * data3[1]
                    elseif length(data3) == size(quad, 4)
                        for (idx4, data4) in enumerate(data3)
                            quad[idx1, idx2, idx3, idx4] = 0.5 * data4
                        end
                    else
                        error("Bad quad size")
                    end
                end
            end
        end
        all(isapprox.(quad, permutedims(quad, [1, 2, 4, 3]); rtol=1e-12)) || error("Asymmetric quad")
        coef_d[2] = quad
    end

    freq, HamiltonianCoefficients{S,M}(coef_d)
end

Base.read(io::IO, ::Type{DenseSystem}) = DenseSystem(read(io, System)...)

Base.read(io::IO, ::Type{DiagonalSystem}) = DiagonalSystem(read(io, System)...)

Base.write(io::IO, sys::System) = JSON.print(io, sys)

function JSON.lower(sys::DenseSystem{S,M}) where {S,M}
    all(sys.freq .== sys.freq[:, 1]) || throw(DomainError(sys.freq, "All surfaces must have the same frequencies."))

    result = Dict()

    result["number of surfaces"] = S
    result["number of modes"] = M

    result["frequencies"] = sys.freq[:, 1]
    if haskey(sys.coef, 0)
        result["energies"] = permutedims(sys.coef[0])
    end
    if haskey(sys.coef, 1)
        result["linear couplings"] = permutedims(sys.coef[1], [3, 2, 1])
    end
    if haskey(sys.coef, 2)
        result["quadratic couplings"] = 2 * permutedims(sys.coef[2], [4, 3, 2, 1])
    end

    result
end

function JSON.lower(sys::DiagonalSystem{S,M}) where {S,M}
    all(sys.freq .== sys.freq[:, 1]) || throw(DomainError(sys.freq, "All surfaces must have the same frequencies."))

    result = Dict()

    result["number of surfaces"] = S
    result["number of modes"] = M

    result["frequencies"] = sys.freq[:, 1]
    if haskey(sys.coef, 0)
        result["energies"] = diag(sys.coef[0])
    end
    if haskey(sys.coef, 1)
        result["linear couplings"] = permutedims(diag(sys.coef[1], 2, 3))
    end
    if haskey(sys.coef, 2)
        result["quadratic couplings"] = 2 * permutedims(diag(sys.coef[2], 3, 4), [3, 2, 1])
    end

    result
end

function show_tensors(io::IO, sys::System{S,M}) where {S,M}
    println(io, "frequencies:")
    for m in 1:M
        println(io, " mode $(m)")
        show_vector(io, sys.freq[m, :])
    end
    for (ord, val) in sys.coef
        if ord == 0
            println(io, "energies:")
        elseif ord == 1
            println(io, "linear couplings:")
        elseif ord == 2
            println(io, "quadratic couplings:")
        else
            println(io, "order $(ord) couplings:")
        end
        for idx in mode_indices(val)
            iszero(val[idx, :, :]) && continue
            if length(idx) == 1
                println(io, " mode $(idx[1])")
            elseif length(idx) >= 2
                println(io, " modes $(join(reverse(Tuple(idx)), ", "))")
            end
            show_matrix(io, val[idx, :, :])
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
