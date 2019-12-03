# Monte Carlo solution for finite temperature with finite difference
# estimators.

"""
    get_sample_fd(sys::System{S,M}, pseudosps::NTuple{N,SamplingParameters{S,M,P}}, sp::SamplingParameters{S_,M,P}, qs::Matrix{Float64})

Compute a sample for `sys` using `pseudosps`, `sp`, and the points `qs`.
"""
function get_sample_fd(sys::System{S,M}, pseudosps::NTuple{N,SamplingParameters{S,M,P}}, sp::SamplingParameters{S_,M,P}, qs::Matrix{Float64}) where {S,S_,M,P,N}
    length(pseudosps) >= 1 || throw(DomainError(length(pseudosps), "At least one set of parameters."))

    # Calculate the numerator and denominator matrices.
    nums = [Matrix{Float64}(I, S, S) for _ in 1:N]
    denom = Matrix{Float64}(I, S_, S_)

    for i1 in 1:P
        # Next index in cyclic path.
        i2 = mod(i1+1, 1:P)

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

    isnothing(sampling_sys) && (sampling_sys = sys_diag_simple)
    sp = SamplingParameters(sampling_sys, beta, P)

    samples_Z = zeros(Float64, num_samples)
    samples_m = zeros(Float64, num_samples)
    samples_p = zeros(Float64, num_samples)

    loop_samples(sp, num_samples, progress_output) do n, qs
        sample = get_sample_fd(sys, (pseudosp, pseudosp_m, pseudosp_p), sp, qs)
        samples_Z[n] = sample[1]
        samples_m[n] = sample[2]
        samples_p[n] = sample[3]
    end

    simple = Analytical(sampling_sys, beta)
    simple_m = Analytical(sampling_sys, beta-dbeta)
    simple_p = Analytical(sampling_sys, beta+dbeta)
    Zrat_m = simple.Z / simple_m.Z
    Zrat_p = simple.Z / simple_p.Z
    normalization = simple.Z

    f_Z(sample_Z) =
        sample_Z * normalization
    f_E(sample_Z, sample_m, sample_p) =
        simple.E +
        1.0/(2dbeta) * (Zrat_m * sample_m - Zrat_p * sample_p) / sample_Z
    f_Cv(sample_Z, sample_m, sample_p) =
        simple.Cv +
        (1.0/dbeta^2 * (Zrat_m * sample_m + Zrat_p * sample_p) / sample_Z -
         2.0/dbeta^2 -
         (E - simple.E)^2) * beta^2

    Z, Z_err = jackknife(f_Z, samples_Z)
    E, E_err = jackknife(f_E, samples_Z, samples_m, samples_p)
    Cv, Cv_err = jackknife(f_Cv, samples_Z, samples_m, samples_p)

    SamplingFiniteDifference(Z, Z_err, E, E_err, Cv, Cv_err)
end
