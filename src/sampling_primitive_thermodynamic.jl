# Monte Carlo solution with primitive thermodynamic estimators.

"""
    get_sample_pt{S,M,P}(sp::SamplingParameters{S,M,P})

Compute a sample using `sp`.
"""
function get_sample_pt{S,M,P}(sp::SamplingParameters{S,M,P})
    # Choose a surface.
    s_ = sample(1:S, sp.weights)

    # Sample coordinates.
    qs = zeros(P, M)
    for m in 1:M
        qs[:, m] = rand(sp.mvns[m, s_])
    end

    # Calculate the numerator, energy, heat capacity, and denominator matrices.
    num = eye(S)
    num_Es = [eye(S) for _ in 1:P]
    num_Cvs = [eye(S) for _ in 1:P]
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

        # Energy and heat capacity estimator matrices.
        EM1 = Symmetric(preIM)
        EM2 = zeros(S, S)
        if i1 == 1
            CM1 = zeros(S, S)
        end
        for s in 1:S
            EM2[s, s] -= sp.sys.energy[s, s] + sp.deltas[s]

            for m in 1:M
                q_disp1 = qs[i1, m] - sp.ds[m, s]
                q_disp2 = qs[i2, m] - sp.ds[m, s]

                EM2[s, s] -= sp.sys.freq[m, s] * 0.5 * sp.Cs[m, s] / sp.Ss[m, s]
                EM2[s, s] -= sp.sys.freq[m, s] * q_disp1 * q_disp2 * sp.Cs[m, s] / sp.Ss[m, s]^2
                EM2[s, s] += sp.sys.freq[m, s] * 0.5 * (q_disp1^2 + q_disp2^2) / sp.Ss[m, s]^2

                if i1 == 1
                    CM1[s, s] += sp.sys.freq[m, s]^2 * 0.5 / sp.Ss[m, s]^2
                    CM1[s, s] += sp.sys.freq[m, s]^2 * 0.5 * q_disp1 * q_disp2 * (sp.Ss[m, s]^2 + sp.Cs[m, s]^2 + 3) / sp.Ss[m, s]^3
                    CM1[s, s] -= sp.sys.freq[m, s]^2 * (q_disp1^2 + q_disp2^2) * sp.Cs[m, s] / sp.Ss[m, s]^3
                end
            end
        end
        EM = EM1 - EM2
        if i1 == 1
            CM = EM1 * EM + CM1 - EM * EM2
        end

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

        IFM = IM * FM
        IEFM = IM * EM * FM

        num *= IM * FM
        for i in 1:P
            if i == i1
                num_Es[i] *= IEFM
                if i1 == 1
                    num_Cvs[i] *= IM * CM * FM
                else
                    num_Cvs[i] *= IEFM
                end
            else
                num_Es[i] *= IFM
                if i1 == 1
                    num_Cvs[i] *= IEFM
                else
                    num_Cvs[i] *= IFM
                end
            end
        end
        denom *= FM
    end

    trace_num = trace(num)
    trace_num_E = mean(trace(num_E) for num_E in num_Es)
    trace_num_Cv = mean(trace(num_Cv) for num_Cv in num_Cvs)
    trace_denom = trace(denom)

    [trace_num/trace_denom, trace_num_E/trace_denom, trace_num_Cv/trace_denom]
end

"""
Monte Carlo solution for an arbitrary system using primitive thermodynamic
estimators.
"""
struct SamplingPrimitiveThermodynamic <: Solution
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
    SamplingPrimitiveThermodynamic{S,M}(sys::System{S,M}, beta::Float64, P::Int, num_samples::Int)

Calculate the solution for `sys` at `beta` with `P` links and `num_samples`
random samples.
"""
function SamplingPrimitiveThermodynamic{S,M}(sys::System{S,M}, beta::Float64, P::Int, num_samples::Int)
    sp = SamplingParameters(sys, beta, P)

    samples = zeros(Float64, 3, num_samples)

    # zero, Inf, NaN
    problems = [false, false, false]

    @showprogress for n in 1:num_samples
        samples[:, n] = get_sample_pt(sp)

        if !problems[1] && samples[1, n] == 0.0
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

    f_E(samples, samples_E, samples_Cv) =
        samples_E ./ samples
    f_Cv(samples, samples_E, samples_Cv) =
        (samples_Cv ./ samples - (samples_E ./ samples).^2) * beta^2

    Z = mean(samples[1, :])
    Z_err = std(samples[1, :]) / sqrt(num_samples)
    E, E_err = jackknife(f_E, [samples[n, :] for n in 1:size(samples, 1)]...)
    Cv, Cv_err = jackknife(f_Cv, [samples[n, :] for n in 1:size(samples, 1)]...)

    SamplingPrimitiveThermodynamic(Z, Z_err, E, E_err, Cv, Cv_err)
end
