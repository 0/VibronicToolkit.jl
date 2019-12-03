# Iterative decomposition of a system into uncoupled surfaces.

"""
The iteration has encountered a degeneracy that it cannot handle.
"""
struct DegeneracyException <: Exception
    n::Int
end

Base.show(io::IO, ex::DegeneracyException) = show(io, "Cannot handle $(ex.n)-fold degeneracy.")

"""
Iterative decomposition of a system into uncoupled surfaces.
"""
struct IterativeDecomposition{S,S_,M}
    "System which was decomposed."
    sys::System{S_,M}

    "Degeneracies at iteration steps."
    degens::Vector{Int}

    "Transformation matrix (S_, S)."
    vs::Matrix{Float64}
    "Coefficients."
    ws::Dict{Int,Array{Float64}}
end

"""
    IterativeDecomposition(sys::System, max_iter::Int; sigma_cutoff::Float64=0.01)

Generate an iterative decomposition of `sys` using no more than `max_iter`
iterations.

Iterations are stopped when the largest singular value drops below the initial
singular value multiplied by `sigma_cutoff`.
"""
function IterativeDecomposition(sys::System{S,M}, max_iter::Int; sigma_cutoff::Float64=0.01) where {S,M}
    all(sys.freq .== sys.freq[:, 1]) || throw(DomainError(sys.freq, "All surfaces must have the same frequencies."))

    # Remaining tensor values.
    xs = Dict{Int,Array{Float64}}(ord => copy(val) for (ord, val) in sys.coef)

    # Generated vectors and coefficients.
    vs = Matrix{Float64}(undef, S, 0)
    ws = Dict{Int,Array{Float64}}(ord => Array{Float64}(undef, repeat([M], ord)..., 0) for (ord, val) in xs)

    degens = Int[]

    S_iter = 0
    first_sigma = nothing
    for _ in 1:max_iter
        # Construct and SVD the matrix of coefficients.
        x = zeros(S, S * sum(M^ord for (ord, val) in xs))
        size(x)[2] > S || error("Too few coefficients")

        col_idx = 0
        for (ord, val) in xs
            for idx in mode_indices(val)
                for s1 in 1:S
                    col_idx += 1
                    for s2 in 1:S
                        x[s2, col_idx] = val[idx, s2, s1]
                    end
                end
            end
        end
        col_idx == size(x)[2] || error("Invalid number of coefficients")
        F = svd(x)
        sigma = F.S[1]

        # Stopping criterion.
        if isnothing(first_sigma)
            first_sigma = sigma
        elseif sigma < sigma_cutoff*first_sigma
            break
        end

        # Determine degeneracy.
        num_degen = 0
        for s in 1:S
            abs(F.S[s] - sigma) / sigma > 1e-15 && break
            num_degen += 1
        end
        push!(degens, num_degen)

        # Compute the vectors.
        pre_vs = Tuple{Vector{Float64},Float64,Vector{Float64},Vector{Float64}}[]
        if num_degen == 1
            v = F.U[:, 1]
            push!(pre_vs, (v, +1, v, v))
        elseif num_degen == 2
            v1 = F.U[:, 1]
            push!(pre_vs, (v1, +1, v1, v1))
            v2 = F.U[:, 2]
            push!(pre_vs, (v2, +1, v2, v2))
            vp = (v1 .+ v2) ./ sqrt(2)
            push!(pre_vs, (vp, +1, v1, v2))
            vm = (v1 .- v2) ./ sqrt(2)
            push!(pre_vs, (vm, -1, v1, v2))
            # Check symmetry.
            for (ord, val) in xs
                for idx in mode_indices(val)
                    v1' * val[idx, :, :] * v2 == v2' * val[idx, :, :] * v1 || error("Asymmetric at order $(ord)")
                end
            end
        else
            throw(DegeneracyException(num_degen))
        end

        # Save vectors and coefficients.
        for (v, s, vL, vR) in pre_vs
            @cat(vs, v)
            for (ord, val) in xs
                w = zeros(repeat([M], ord)...)
                for idx in mode_indices(val)
                    w[idx] = s * vL' * val[idx, :, :] * vR
                end
                ws[ord] = cat(ws[ord], w; dims=length(size(ws[ord])))
            end
            S_iter += 1
        end

        # Update arrays.
        for n in 1:length(pre_vs)
            v = pre_vs[n][1]
            S_idx = size(vs, 2) - length(pre_vs) + n
            for (ord, val) in xs
                for idx in mode_indices(val)
                    val[idx, :, :] .-= ws[ord][idx, S_idx] * v * v'
                end
            end
        end
    end

    # Check that we recover the original tensors.
    for (ord, val) in xs
        for idx in mode_indices(val)
            recovered = val[idx, :, :] + sum([ws[ord][idx, n]*vs[:, n]*vs[:, n]' for n in 1:S_iter])
            err = maximum(abs.(sys.coef[ord][idx, :, :] - recovered))
            err < 1e-14 || error("Maximum deviation at order $(ord): $(err)")
        end
    end

    IterativeDecomposition{S_iter,S,M}(sys, degens, vs, ws)
end

"""
    DiagonalSystem(decomp::IterativeDecomposition)

Generate a `DiagonalSystem` from `decomp`.
"""
function DiagonalSystem(decomp::IterativeDecomposition{S,S_,M}) where {S,S_,M}
    freq = zeros(M, S)
    for s in 1:S
        for m in 1:M
            freq[m, s] = decomp.sys.freq[m, 1]
        end
    end

    coef = Dict{Int,Array{Float64}}(ord => zeros(size(val)..., S) for (ord, val) in decomp.ws)
    for (ord, val) in decomp.ws
        for s in 1:S
            for idx in mode_indices(val; dims=1)
                coef[ord][idx, s, s] = val[idx, s]
            end
        end
    end
    DiagonalSystem(freq, coef)
end
