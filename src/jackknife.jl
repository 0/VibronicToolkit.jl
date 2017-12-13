# Jackknife resampling.

"""
    jackknife(f, xss...)

Perform jackknife resampling for the function `f` with the data `xss`.

This uses the formulas from Peter Young's notes
(https://arxiv.org/abs/1210.3781, https://arxiv.org/pdf/1210.3781v3.pdf).
"""
function jackknife(f, xss...)
    # At least one data set.
    length(xss) >= 1 || throw(DomainError())

    N = length(xss[1])
    # Uniform length.
    all(length(xs) .== N for xs in xss) || throw(DomainError())

    # Regular averages.
    xss_mean = Float64[]
    for xs in xss
        push!(xss_mean, mean(xs))
    end

    # "All-but-one" averages.
    xss_J = Array{Float64,1}[]
    for xs in xss
        push!(xss_J, (sum(xs) - xs) / (N - 1))
    end

    f_plain = f(xss_mean...)
    f_Js = f(xss_J...)

    # Unbiased mean.
    f_J = N * f_plain - (N - 1) * mean(f_Js)
    # Standard error.
    f_err = sqrt(N - 1) * std(f_Js; corrected=false)

    f_J, f_err
end
