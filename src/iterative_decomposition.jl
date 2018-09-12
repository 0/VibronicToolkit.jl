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

    "Coefficients for the energy (S)."
    ws_energy::Vector{Float64}
    "Coefficients for the linear coupling (M, S)."
    ws_lin::Matrix{Float64}
    "Coefficients for the quadratic coupling (M, M, S)."
    ws_quad::Array{Float64,3}
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
    x_energy = copy(sys.energy)
    x_lin = copy(sys.lin)
    x_quad = copy(sys.quad)

    # Generated vectors and coefficients.
    vs = Matrix{Float64}(undef, S, 0)
    ws_energy = Vector{Float64}(undef, 0)
    ws_lin = Matrix{Float64}(undef, M, 0)
    ws_quad = Array{Float64,3}(undef, M, M, 0)

    degens = Int[]

    first_sigma = nothing
    for _ in 1:max_iter
        # Construct and SVD the matrix of coefficients.
        x = zeros(S, S + S*M + S*M*M)
        for s1 in 1:S
            row = s1
            for s2 in 1:S
                col0 = s2
                x[row, col0] = x_energy[s1, s2]
                for m1 in 1:M
                    col1 = S + (s2-1)*M + m1
                    x[row, col1] = x_lin[m1, s1, s2]
                    for m2 in 1:M
                        col2 = S + S*M + (s2-1)*M*M + (m1-1)*M + m2
                        x[row, col2] = x_quad[m1, m2, s1, s2]
                    end
                end
            end
        end
        F = svd(x)
        sigma = F.S[1]

        # Stopping criterion.
        if first_sigma === nothing
            first_sigma = sigma
        elseif sigma < sigma_cutoff*first_sigma
            break
        end

        # Determine degeneracy.
        num_degen = 0
        for s in 1:S
            F.S[s] != sigma && break
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
            v1' * x_energy * v2 == v2' * x_energy * v1 || error("Asymmetric x_energy")
            for m1 in 1:M
                v1' * x_lin[m1, :, :] * v2 == v2' * x_lin[m1, :, :] * v1 || error("Asymmetric x_lin")
                for m2 in 1:M
                    v1' * x_quad[m2, m1, :, :] * v2 == v2' * x_quad[m2, m1, :, :] * v1 || error("Asymmetric x_quad")
                end
            end
        else
            throw(DegeneracyException(num_degen))
        end

        # Save vectors and coefficients.
        for (v, s, vL, vR) in pre_vs
            w_energy = s * vL' * x_energy * vR
            w_lin = zeros(M)
            w_quad = zeros(M, M)
            for m1 in 1:M
                w_lin[m1] = s * vL' * x_lin[m1, :, :] * vR
                for m2 in 1:M
                    w_quad[m2, m1] = s * vL' * x_quad[m2, m1, :, :] * vR
                end
            end
            @cat(vs, v)
            @cat(ws_energy, w_energy)
            @cat(ws_lin, w_lin)
            @cat(ws_quad, w_quad)
        end

        # Update arrays.
        for n in 1:length(pre_vs)
            v = pre_vs[n][1]
            idx = size(vs, 2) - length(pre_vs) + n
            x_energy .-= ws_energy[idx] * v * v'
            for m1 in 1:M
                x_lin[m1, :, :] .-= ws_lin[m1, idx] * v * v'
                for m2 in 1:M
                    x_quad[m2, m1, :, :] .-= ws_quad[m2, m1, idx] * v * v'
                end
            end
        end
    end

    # Number of generated surfaces.
    S_iter = length(ws_energy)

    # Check that we recover the original tensors.
    recovered_energy = x_energy + sum([ws_energy[n]*vs[:, n]*vs[:, n]' for n in 1:S_iter])
    err = maximum(abs.((sys.energy - recovered_energy)))
    err < 1e-14 || error("Maximum deviation: $(err)")
    for m1 in 1:M
        recovered_lin = x_lin[m1, :, :] + sum([ws_lin[m1, n]*vs[:, n]*vs[:, n]' for n in 1:S_iter])
        err = maximum(abs.((sys.lin[m1, :, :] - recovered_lin)))
        err < 1e-14 || error("Maximum deviation: $(err)")
        for m2 in 1:M
            recovered_quad = x_quad[m2, m1, :, :] + sum([ws_quad[m2, m1, n]*vs[:, n]*vs[:, n]' for n in 1:S_iter])
            err = maximum(abs.((sys.quad[m2, m1, :, :] - recovered_quad)))
            err < 1e-14 || error("Maximum deviation: $(err)")
        end
    end

    IterativeDecomposition{S_iter,S,M}(sys, degens, vs, ws_energy, ws_lin, ws_quad)
end

"""
    DiagonalSystem(decomp::IterativeDecomposition)

Generate a `DiagonalSystem` from `decomp`.
"""
function DiagonalSystem(decomp::IterativeDecomposition{S,S_,M}) where {S,S_,M}
    energy = zeros(S, S)
    freq = zeros(M, S)
    lin = zeros(M, S, S)
    quad = zeros(M, M, S, S)
    for s in 1:S
        energy[s, s] = decomp.ws_energy[s]
        for m1 in 1:M
            freq[m1, s] = decomp.sys.freq[m1, 1]
            lin[m1, s, s] = decomp.ws_lin[m1, s]
            for m2 in 1:M
                quad[m1, m2, s, s] = decomp.ws_quad[m1, m2, s]
            end
        end
    end
    DiagonalSystem(energy, freq, lin, quad)
end
