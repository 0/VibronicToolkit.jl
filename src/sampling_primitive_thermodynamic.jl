# Monte Carlo solution with primitive thermodynamic estimators.

"""
    get_sample_pt(sys::System{S,M}, sp::SamplingParameters{S,M,P})

Compute a sample for `sys` using `sp`.
"""
function get_sample_pt(sys::System{S,M}, sp::SamplingParameters{S,M,P}) where {S,M,P}
    # Choose a surface.
    s_ = sample(1:S, sp.weights)

    # Sample coordinates.
    qs = zeros(P, M)
    for m in 1:M
        qs[:, m] .= rand(sp.mvns[m, s_])
    end

    # Calculate the numerator, energy, heat capacity, and denominator matrices.
    num = Matrix{Float64}(I, S, S)
    num_Es = [Matrix{Float64}(I, S, S) for _ in 1:P]
    num_Cvs = [Matrix{Float64}(I, S, S) for _ in 1:P]
    denom = Matrix{Float64}(I, S, S)

    for i1 in 1:P
        # Next index in cyclic path.
        i2 = i1%P+1

        FM, _ = sampling_matrix_free_particle(sp, qs[i1, :], qs[i2, :])
        EM1, IM = sampling_matrix_interaction(sys, sp, qs[i1, :])
        EM2 = sampling_matrix_energy(sp, qs[i1, :], qs[i2, :])
        EM = EM1 - EM2
        if i1 == 1
            CM1 = sampling_matrix_heat_capacity(sp, qs[i1, :], qs[i2, :])
            CM = EM1 * EM + CM1 - EM * EM2
        end

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

    trace_num = tr(num)
    trace_num_E = mean(tr(num_E) for num_E in num_Es)
    trace_num_Cv = mean(tr(num_Cv) for num_Cv in num_Cvs)
    trace_denom = tr(denom)

    [trace_num/trace_denom, trace_num_E/trace_denom, trace_num_Cv/trace_denom]
end

"""
Monte Carlo solution for an arbitrary system using primitive thermodynamic
estimators.
"""
struct SamplingPrimitiveThermodynamic <: Sampling
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
    SamplingPrimitiveThermodynamic(sys::System{S,M}, beta::Float64, P::Int, num_samples::Int; progress_output::IO=stderr)

Calculate the solution for `sys` at `beta` with `P` links and `num_samples`
random samples.

The progress meter is written to `progress_output`.
"""
function SamplingPrimitiveThermodynamic(sys::System{S,M}, beta::Float64, P::Int, num_samples::Int; progress_output::IO=stderr) where {S,M}
    sys_diag = diag(sys)
    sys_diag_simple = simplify(sys_diag)

    sp = SamplingParameters(sys_diag_simple, beta, P)

    samples = zeros(Float64, 3, num_samples)

    # zero, Inf, NaN
    problems = [false, false, false]

    meter = Progress(num_samples, output=progress_output)
    for n in ProgressWrapper(1:num_samples, meter)
        samples[:, n] .= get_sample_pt(sys, sp)

        if !problems[1] && samples[1, n] == 0.0
            @warn "\azero!"
            problems[1] = true
        end

        if !problems[2] && any(samples[:, n] .== Inf)
            @warn "\aInf!"
            problems[2] = true
        end

        if !problems[3] && any(samples[:, n] .!= samples[:, n])
            @warn "\aNaN!"
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
