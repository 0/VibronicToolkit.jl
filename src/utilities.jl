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
