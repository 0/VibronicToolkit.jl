# Ground state potential energy surface.

"""
    ground_pes(sys::System{S,2}, extent::NTuple{4,Float64}; lengths::NTuple{2,Float64}=(501, 501))

Compute the ground state PES for the 2-mode system `sys`.

The extent for each mode is given in `extent`, and the number of points is
given in `lengths`.
"""
function ground_pes(sys::System{S,2}, extent::NTuple{4,Float64}; lengths::NTuple{2,Int}=(501, 501)) where {S}
    M = 2
    q1_min, q1_max, q2_min, q2_max = extent
    q1_length, q2_length = lengths
    pes = zeros(q2_length, q1_length)
    for (i, q1) in enumerate(range(q1_min; stop=q1_max, length=q1_length))
        for (j, q2) in enumerate(range(q2_min; stop=q2_max, length=q2_length))
            # Lowest eigenvalue of potential matrix.
            pes[j, i] = eigvals(potential(sys, [q1, q2]))[1]
        end
    end
    pes
end

"""
    path_mean_std(sys::DiagonalSystem, beta::Float64, P::Int, s::Int, m::Int)

Compute the standard deviation of the mean of the path for surface `s` and mode
`m` in `sp`.
"""
function path_mean_std(sys::DiagonalSystem{S,M}, beta::Float64, P::Int, s::Int, m::Int) where {S,M}
    1 <= s <= S || throw(DomainError(s, "Invalid surface index."))
    1 <= m <= M || throw(DomainError(m, "Invalid mode index."))
    sp = SamplingParameters(sys, beta, P)
    centroid_std = sqrt(1 / eigvals(sp.precs[m, s])[1])
    centroid_std / sqrt(P)
end
