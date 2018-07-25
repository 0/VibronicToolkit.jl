# Monte Carlo solution with finite difference estimators.

"""
    get_sample_fd{S,M,P}(sps::SamplingParameters{S,M,P}...)

Compute a sample using `sps`.
"""
function get_sample_fd{S,M,P}(sps::SamplingParameters{S,M,P}...)
    # At least one set of parameters.
    length(sps) >= 1 || throw(DomainError())

    # Choose a surface.
    s_ = sample(1:S, sps[1].weights)

    # Sample coordinates.
    qs = zeros(P, M)
    for m in 1:M
        qs[:, m] = rand(sps[1].mvns[m, s_])
    end

    traces_num = Float64[]
    traces_denom = Float64[]

    for sp in sps
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

            num *= IM * FM
            denom *= FM
        end

        push!(traces_num, trace(num))
        push!(traces_denom, trace(denom))
    end

    [traces_num[1]/traces_denom[1],
     traces_num[2]/traces_denom[1],
     traces_num[3]/traces_denom[1]]
end

"""
Monte Carlo solution for an arbitrary system using finite difference
estimators.
"""
struct SamplingFiniteDifference <: Sampling
    "Partition function."
    Z::Float64
    "Partition function standard error."
    Z_err::Float64
    "Energy."
    E::Float64
    "Energy standard error."
    E_err::Float64
    "Heat capacity."
    Cv::Float64
    "Heat capacity standard error."
    Cv_err::Float64
end

"""
    SamplingFiniteDifference{S,M}(sys::System{S,M}, beta::Float64, dbeta::Float64, P::Int, num_samples::Int)

Calculate the solution for `sys` at `beta` with `P` links and `num_samples`
random samples, using finite difference step `dbeta`.
"""
function SamplingFiniteDifference{S,M}(sys::System{S,M}, beta::Float64, dbeta::Float64, P::Int, num_samples::Int)
    simple = Analytical(simplify(sys), beta)
    simple_m = Analytical(simplify(sys), beta-dbeta)
    simple_p = Analytical(simplify(sys), beta+dbeta)

    Zrat_m = simple.Z / simple_m.Z
    Zrat_p = simple.Z / simple_p.Z

    sp = SamplingParameters(sys, beta, P)
    sp_m = SamplingParameters(sys, beta-dbeta, P)
    sp_p = SamplingParameters(sys, beta+dbeta, P)

    samples = zeros(Float64, 3, num_samples)

    # zero, Inf, NaN
    problems = [false, false, false]

    @showprogress for n in 1:num_samples
        samples[:, n] = get_sample_fd(sp, sp_m, sp_p)

        if !problems[1] && any(samples[:, n] .== 0.0)
            warn("\azero!")
            problems[1] = true
        end

        if !problems[2] && any(samples[:, n] .== Inf)
            warn("\aInf!")
            problems[2] = true
        end

        if !problems[3] && any(samples[:, n] .!= samples[:, n])
            warn("\aNaN!")
            problems[3] = true
        end
    end

    f_E(samples, samples_m, samples_p) =
        simple.E +
        1.0/(2dbeta) * (Zrat_m * samples_m .- Zrat_p * samples_p) ./ samples
    f_Cv(samples, samples_m, samples_p) =
        simple.Cv +
        (1.0/dbeta^2 * (Zrat_m * samples_m .+ Zrat_p * samples_p) ./ samples -
         2.0/dbeta^2 -
         (E - simple.E)^2) * beta^2

    Z = mean(samples[1, :])
    Z_err = std(samples[1, :]) / sqrt(num_samples)
    E, E_err = jackknife(f_E, [samples[n, :] for n in 1:size(samples, 1)]...)
    Cv, Cv_err = jackknife(f_Cv, [samples[n, :] for n in 1:size(samples, 1)]...)

    SamplingFiniteDifference(Z, Z_err, E, E_err, Cv, Cv_err)
end
