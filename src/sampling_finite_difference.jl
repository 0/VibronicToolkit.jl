# Monte Carlo solution for finite temperature with finite difference
# estimators.

"""
    get_sample_fd(sys::System{S,M}, pseudosps::NTuple{N,SamplingParameters{S,M,P}}, sp::SamplingParameters{S_,M,P})

Compute a sample for `sys` using `pseudosps` and `sp`.
"""
function get_sample_fd(sys::System{S,M}, pseudosps::NTuple{N,SamplingParameters{S,M,P}}, sp::SamplingParameters{S_,M,P}) where {S,S_,M,P,N}
    length(pseudosps) >= 1 || throw(DomainError(length(pseudosps), "At least one set of parameters."))

    # Choose a surface.
    s_ = sample(1:S_, sp.weights)

    # Sample coordinates.
    qs = zeros(P, M)
    for m in 1:M
        qs[:, m] .= rand(sp.mvns[m, s_])
    end

    # Calculate the numerator and denominator matrices.
    nums = [Matrix{Float64}(I, S, S) for _ in 1:N]
    denom = Matrix{Float64}(I, S_, S_)

    for i1 in 1:P
        # Next index in cyclic path.
        i2 = i1%P+1

        scaling = nothing

        # Numerator.
        for (j, pseudosp) in enumerate(pseudosps)
            FM_num, scaling = sampling_matrix_free_particle(pseudosp, qs[i1, :], qs[i2, :]; scaling=scaling)
            _, IM = sampling_matrix_interaction(sys, pseudosp, qs[i1, :])
            nums[j] *= IM * FM_num
        end

        # Denominator.
        FM_denom, _ = sampling_matrix_free_particle(sp, qs[i1, :], qs[i2, :]; scaling=scaling)
        denom *= FM_denom
    end

    traces_num = [tr(num) for num in nums]
    trace_denom = tr(denom)

    [trace_num/trace_denom for trace_num in traces_num]
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
    SamplingFiniteDifference(sys::System, beta::Float64, dbeta::Float64, P::Int, num_samples::Int; sampling_sys::Maybe{System}=nothing, progress_output::IO=stderr)

Calculate the solution for `sys` at `beta` using `P` links, `num_samples`
random samples, and finite difference step `dbeta`.

If `sampling_sys` is provided, it is used for sampling. Otherwise, sampling
defaults to the simplified diagonal subsystem of `sys`.

The progress meter is written to `progress_output`.
"""
function SamplingFiniteDifference(sys::System, beta::Float64, dbeta::Float64, P::Int, num_samples::Int; sampling_sys::Maybe{System}=nothing, progress_output::IO=stderr)
    sys_diag = diag(sys)
    sys_diag_simple = simplify(sys_diag)
    pseudosp = SamplingParameters(sys_diag_simple, beta, P)
    pseudosp_m = SamplingParameters(sys_diag_simple, beta-dbeta, P)
    pseudosp_p = SamplingParameters(sys_diag_simple, beta+dbeta, P)

    sampling_sys === nothing && (sampling_sys = sys_diag_simple)
    sp = SamplingParameters(sampling_sys, beta, P)

    samples = zeros(Float64, 3, num_samples)
    meter = Progress(num_samples, output=progress_output)
    for n in ProgressWrapper(1:num_samples, meter)
        samples[:, n] .= get_sample_fd(sys, (pseudosp, pseudosp_m, pseudosp_p), sp)
    end

    simple = Analytical(sampling_sys, beta)
    simple_m = Analytical(sampling_sys, beta-dbeta)
    simple_p = Analytical(sampling_sys, beta+dbeta)
    Zrat_m = simple.Z / simple_m.Z
    Zrat_p = simple.Z / simple_p.Z
    normalization = simple.Z

    f_E(sample, sample_m, sample_p) =
        simple.E +
        1.0/(2dbeta) * (Zrat_m * sample_m - Zrat_p * sample_p) / sample
    f_Cv(sample, sample_m, sample_p) =
        simple.Cv +
        (1.0/dbeta^2 * (Zrat_m * sample_m + Zrat_p * sample_p) / sample -
         2.0/dbeta^2 -
         (E - simple.E)^2) * beta^2

    Z = mean(samples[1, :]) * normalization
    Z_err = std(samples[1, :]) / sqrt(num_samples) * normalization
    E, E_err = jackknife(f_E, [samples[n, :] for n in 1:size(samples, 1)]...)
    Cv, Cv_err = jackknife(f_Cv, [samples[n, :] for n in 1:size(samples, 1)]...)

    SamplingFiniteDifference(Z, Z_err, E, E_err, Cv, Cv_err)
end
