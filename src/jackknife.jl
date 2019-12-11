# Jackknife resampling.

"""
    jackknife(f, xss...; plain_delta::Real=0)

Perform jackknife resampling for the function `f` with the data `xss`.

The contribution of the regular average to the unbiased mean is modulated by
`plain_delta`.

This uses the formulas from Peter Young's notes
(https://arxiv.org/abs/1210.3781, https://arxiv.org/pdf/1210.3781v3.pdf).
"""
function jackknife(f, xss...; plain_delta::Real=0)
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
    f_J = (N + plain_delta) * f_plain - (N - 1) * mean(f_Js)
    # Standard error.
    f_err = sqrt(N - 1) * std(f_Js; corrected=false)

    f_J, f_err
end

"""
    jackknife(f, ws::Weights, labels::Vector{Int}, xss...)

Perform jackknife resampling for the function `f` with the data `xss` using
`ws` and the corresponding `labels`.
"""
function jackknife(f, ws::Weights, labels::Vector{Int}, xss...)
    N = length(labels)
    lengths = [length(xs) for xs in xss]
    all(lengths .== N) || throw(DomainError((N, lengths), "Uniform length, matching labels."))

    S = length(ws)
    nonempty = [s in labels for s in 1:S]
    any(nonempty) || throw(DomainError(nonempty, "At least one non-empty surface."))
    isempty(setdiff(labels, 1:S)) || throw(DomainError(labels, "Only valid labels."))

    protos = [zero(xs[1]) for xs in xss]
    xsss = [[xs[labels .== s] for xs in xss] for s in 1:S]

    pre_means = Array{Any}(undef, length(xss), S)
    for s in 1:S
        for i in 1:length(xss)
            if nonempty[s]
                pre_means[i, s] = ws[s] * mean(xsss[s][i]) / ws.sum
            else
                pre_means[i, s] = protos[i]
            end
        end
    end

    means = zeros(Float64, S)
    errs = zeros(Float64, S)
    for s in 1:S
        nonempty[s] || continue
        function g(ys...)
            xs = copy(protos)
            for s_ in 1:S
                if s_ == s
                    xs .+= ws[s_] .* ys ./ ws.sum
                else
                    xs .+= pre_means[:, s_]
                end
            end
            f(xs...)
        end
        means[s], errs[s] = jackknife(g, xsss[s]...; plain_delta=-1)
    end

    combined_pre_means = copy(protos)
    for s in 1:S
        combined_pre_means .+= pre_means[:, s]
    end

    # Unbiased mean.
    f_J = f(combined_pre_means...) + sum(means)
    # Standard error.
    f_err = sqrt(sum(errs.^2))

    f_J, f_err
end
