# Diagonal density matrix.

"""
    diagonal_density(sys::System{S,2}, beta::Float64, basis_size::Int, extent::NTuple{4,Float64}; lengths::NTuple{2,Int}=(101, 101), progress_output::IO=stderr)

Compute the diagonal (in position) density matrix for the 2-mode system `sys`
at `beta` using `basis_size` basis functions.

The extent for each mode is given in `extent`, and the number of points is
given in `lengths`.

The progress meter is written to `progress_output`.
"""
function diagonal_density(sys::System{S,2}, beta::Float64, basis_size::Int, extent::NTuple{4,Float64}; lengths::NTuple{2,Int}=(101, 101), progress_output::IO=stderr) where {S}
    basis = Basis(sys, basis_size)
    h0, V = operators(basis, sys)
    basis_vectors = vectors(basis)

    F = eigen(Symmetric(h0 + V))
    Es = F.values
    Vs = F.vectors

    Z = sum(exp.(-beta * Es))

    M = 2
    q1_min, q1_max, q2_min, q2_max = extent
    q1_length, q2_length = lengths
    density = zeros(q2_length, q1_length, S, S)
    wfs = Array{Float64}(undef, basis.dim, S)
    meter = Progress(q1_length, output=progress_output)
    for (i, q1) in ProgressWrapper(enumerate(range(q1_min; stop=q1_max, length=q1_length)), meter)
        for (j, q2) in enumerate(range(q2_min; stop=q2_max, length=q2_length))
            fill!(wfs, 0.0)
            for idx in 1:basis.dim
                s = basis_vectors[M+1, idx]
                wfs[:, s] .+= ho_wf(basis_vectors[1:M, idx], [q1, q2]) .* Vs[idx, :]
            end
            for s1 in 1:S
                for s2 in 1:S
                    density[j, i, s2, s1] = sum(exp.(-beta * Es) .* wfs[:, s1] .* wfs[:, s2]) / Z
                end
            end
        end
    end
    density
end
