# Assorted utility functions.

"""
    diag(A::Array, dim1::Int, dim2::Int)

Collapse `dim1` and `dim2` of `A` into a single dimension containing only the
elements at equal indices.
"""
function diag(A::Array, dim1::Int, dim2::Int)
    s = [size(A)...]
    s[dim1] == s[dim2] || throw(DomainError(s, "Dimensions must have the same size."))

    deleteat!(s, dim2)

    result = similar(A, s...)
    for idx1 in CartesianIndices(A)
        idx1[dim1] == idx1[dim2] || continue
        idx2 = [idx1.I...]
        deleteat!(idx2, dim2)
        result[idx2...] = A[idx1]
    end

    result
end

"""
    hermite(n::Int, q::Float64)

Compute the value `H_n(q)/sqrt(sqrt(pi) 2^n n!)` of the normalized `n`th
Hermite polynomial at `q`.
"""
function hermite(n::Int, q::Float64)
    h1 = 1.0 / sqrt(sqrt(pi))
    n == 0 && return h1
    h2 = 2q / sqrt(2*sqrt(pi))
    for m in 1:(n-1)
        h1, h2 = h2, (sqrt(2) * q * h2 - sqrt(m) * h1) / sqrt(m+1)
    end
    h2
end

"""
    ho_wf(n::Int, q::Float64)

Compute the value of the `n`th harmonic oscillator wavefunction at `q`.
"""
ho_wf(n::Int, q::Float64) = exp(-q^2/2) * hermite(n, q)

"""
    ho_wf(ns::Vector{Int}, qs::Vector{Float64})

Compute the value of the multi-mode harmonic oscillator wavefunction.
"""
ho_wf(ns::Vector{Int}, qs::Vector{Float64}) = prod(x -> ho_wf(x...), zip(ns, qs))

"""
    gauss_cdf_inv(phi::Float64)

Invert the standard Gaussian cdf `phi`.
"""
gauss_cdf_inv(phi::Float64) = -sqrt(2) * erfcinv(2 * phi)

"""
    @cat A B

Append `B` along the last dimension of `A`, modifying the binding of `A`.
"""
macro cat(A::Symbol, B)
    quote
        $(esc(A)) = cat($(esc(A)), $(esc(B)); dims=length(size($(esc(A)))))
    end
end

Maybe{T} = Union{T,Nothing} where {T}

function show_value(io::IO, x::Float64)
    if x == 0.0
        print(io, "    .          ")
    else
        @printf(io, "% 15.10f", x)
    end
    nothing
end

function show_vector(io::IO, xs::Vector{Float64})
    for i in 1:length(xs)
        show_value(io, xs[i])
    end
    println(io)
    nothing
end

function show_matrix(io::IO, xs::Matrix{Float64})
    for i in 1:size(xs, 2)
        show_vector(io, xs[:, i])
    end
    nothing
end
