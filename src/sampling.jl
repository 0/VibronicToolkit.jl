# Monte Carlo solution for arbitrary systems.

"""
Sampling parameters for a particular system.
"""
struct SamplingParameters{S,M,P}
    sys::System{S,M}

    "Imaginary time step."
    tau::Float64

    "Surface weights."
    weights::WeightVec
    "Multivariate normal distributions."
    mvns::Matrix{MvNormal}

    "Energy offsets due to linear terms."
    deltas::Vector{Float64}
    "Position offsets due to linear terms."
    ds::Matrix{Float64}

    "Precomputed cosh factors."
    Cs::Matrix{Float64}
    "Precomputed sinh factors."
    Ss::Matrix{Float64}
    "Precomputed products of sinh factors."
    S_prods::Vector{Float64}
end

"""
    SamplingParameters{S,M}(sys::System{S,M}, beta::Float64, P::Int)

Generate sampling parameters for `sys` at `beta` with `P` links.
"""
function SamplingParameters{S,M}(sys::System{S,M}, beta::Float64, P::Int)
    tau = beta / P

    Zas = ones(S)
    deltas = zeros(S)
    ds = zeros(M, S)
    for s in 1:S
        for m in 1:M
            Zas[s] *= 1.0 / (2*sinh(0.5 * beta * sys.freq[m, s]))
            deltas[s] += -0.5 * sys.lin[m, s, s].^2 / sys.freq[m, s]
            ds[m, s] = -sys.lin[m, s, s] / sys.freq[m, s]
        end
    end

    preweights = exp.(-beta * (diag(sys.energy) + deltas))
    weights = WeightVec(preweights .* Zas)

    Cs = cosh.(sys.freq * tau)
    Ss = sinh.(sys.freq * tau)
    S_prods = Float64[prod(Ss[:, s]) for s in 1:S]

    mvns = Matrix{MvNormal}(M, S)
    for s in 1:S
        for m in 1:M
            mean = ds[m, s] * ones(P)

            # Precision matrix.
            prec = 2 * Cs[m, s] * eye(P) - diagm(ones(P-1), -1) - diagm(ones(P-1), 1)
            prec[1, end] = -1.0
            prec[end, 1] = -1.0

            # Convariance matrix.
            cov = inv(prec / Ss[m, s])
            maximum(abs.(cov' - cov)) < 1e-13 || warn("Asymmetric cov: $(maximum(abs.(cov' - cov)))")
            # Force it to be symmetric.
            cov = (cov + cov') / 2

            mvns[m, s] = MvNormal(mean, cov)
        end
    end

    SamplingParameters{S,M,P}(sys, tau, weights, mvns, deltas, ds, Cs, Ss, S_prods)
end

"""
    get_sample{S,M,P}(sp::SamplingParameters{S,M,P})

Compute a sample using `sp`.
"""
function get_sample{S,M,P}(sp::SamplingParameters{S,M,P})
    # Choose a surface.
    s_ = sample(1:S, sp.weights)

    # Sample coordinates.
    qs = zeros(P, M)
    for m in 1:M
        qs[:, m] = rand(sp.mvns[m, s_])
    end

    # Calculate the numerator and denominator matrices.
    num = eye(S)
    denom = eye(S)

    for i1 in 1:P
        # Next index in cyclic path.
        i2 = i1%P+1

        # Interaction matrix.
        preIM = zeros(S, S)
        for s1 in 1:S
            # Only build the upper triangle.
            for s2 in 1:s1
                if s1 != s2
                    preIM[s2, s1] += sp.sys.energy[s2, s1]

                    for m in 1:M
                        preIM[s2, s1] += sp.sys.lin[m, s2, s1] * qs[i1, m]
                    end
                end

                for m1 in 1:M
                    for m2 in 1:M
                        preIM[s2, s1] += 0.5 * sp.sys.quad[m2, m1, s2, s1] * qs[i1, m2] * qs[i1, m1]
                    end
                end
            end
        end
        IM = expm(Symmetric(-sp.tau * preIM))

        # Free particle matrix.
        FM = zeros(S, S)
        for s in 1:S
            x = -sp.tau * (sp.sys.energy[s, s] + sp.deltas[s])

            for m in 1:M
                q_disp1 = qs[i1, m] - sp.ds[m, s]
                q_disp2 = qs[i2, m] - sp.ds[m, s]

                x += -0.5 * ((q_disp1^2 + q_disp2^2) * sp.Cs[m, s] - 2 * q_disp1 * q_disp2) / sp.Ss[m, s]
            end

            FM[s, s] += exp(x) / sqrt(sp.S_prods[s])
        end
        FM /= maximum(FM)

        num *= IM * FM
        denom *= FM
    end

    trace(num) / trace(denom)
end

"""
Monte Carlo solution for an arbitrary system.
"""
struct Sampling <: Solution
    "Partition function."
    Z::Float64
    "Partition function standard error."
    Z_err::Float64
end

"""
    Sampling{S,M}(sys::System{S,M}, beta::Float64, P::Int, num_samples::Int)

Calculate the solution for `sys` at `beta` with `P` links and `num_samples`
random samples.
"""
function Sampling{S,M}(sys::System{S,M}, beta::Float64, P::Int, num_samples::Int)
    sp = SamplingParameters(sys, beta, P)

    samples = zeros(num_samples)
    # zero, Inf, NaN
    problems = [false, false, false]

    @showprogress for n in 1:num_samples
        samples[n] = get_sample(sp)

        if !problems[1] && samples[n] == 0.0
            warn("\azero!")
            problems[1] = true
        end

        if !problems[2] && samples[n] == Inf
            warn("\aInf!")
            problems[2] = true
        end

        if !problems[3] && samples[n] != samples[n]
            warn("\aNaN!")
            problems[3] = true
        end
    end

    Z = mean(samples)
    Z_err = std(samples) / sqrt(length(samples))

    Sampling(Z, Z_err)
end
