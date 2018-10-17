# Jackknife resampling.

"""
    jackknife(f, xss...)

Perform jackknife resampling for the function `f` with the data `xss`.

This uses the formulas from Peter Young's notes
(https://arxiv.org/abs/1210.3781, https://arxiv.org/pdf/1210.3781v3.pdf).
"""
function jackknife(f, xss...)
    length(xss) >= 1 || throw(DomainError(length(xss), "At least one data set."))

    N = length(xss[1])
    lengths = [length(xs) for xs in xss]
    all(lengths .== N) || throw(DomainError(lengths, "Uniform length."))

    # Regular averages.
    xss_mean = [mean(xs) for xs in xss]
    f_plain = f(xss_mean...)

    # "All-but-one" averages.
    sums = [sum(xs) for xs in xss]
    f_Js = Vector{Float64}(undef, N)
    for n in 1:N
        xss_J = [(s - xs[n]) / (N - 1) for (s, xs) in zip(sums, xss)]
        f_Js[n] = f(xss_J...)
    end

    # Unbiased mean.
    f_J = N * f_plain - (N - 1) * mean(f_Js)
    # Standard error.
    f_err = sqrt(N - 1) * std(f_Js; corrected=false)

    f_J, f_err
end
