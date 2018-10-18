# Monte Carlo solution for PIGS.

"""
    get_sample_pigs(sys::System{S,M}, pseudosp::PigsSamplingParameters{S,M,P}, sp::PigsSamplingParameters{S_,M,P})

Compute a sample for `sys` using `pseudosp` and `sp`.
"""
function get_sample_pigs(sys::System{S,M}, pseudosp::PigsSamplingParameters{S,M,P}, sp::PigsSamplingParameters{S_,M,P}) where {S,S_,M,P}
    # Choose a surface.
    s_ = sample(1:S_, sp.weights)

    # Sample coordinates.
    qs = zeros(P+1, M)
    for m in 1:M
        qs[:, m] .= rand(sp.mvns[m, s_])
    end

    # Calculate the numerator and denominator matrices.
    num = Matrix{Float64}(I, S, S)
    denom = Matrix{Float64}(I, S_, S_)

    for i1 in 1:P
        # Next index in open path.
        i2 = i1+1

        # Numerator.
        FM_num, scaling = sampling_matrix_free_particle(pseudosp, qs[i1, :], qs[i2, :])
        _, IM_L = sampling_matrix_interaction(sys, pseudosp, qs[i1, :])
        _, IM_R = sampling_matrix_interaction(sys, pseudosp, qs[i2, :])
        num *= sqrt(IM_L) * FM_num * sqrt(IM_R)

        # Denominator.
        FM_denom, _ = sampling_matrix_free_particle(sp, qs[i1, :], qs[i2, :]; scaling=scaling)
        denom *= FM_denom
    end

    # Average the energy at the two ends.
    num_E = 0.5 * (potential(sys, qs[1, :]) * num + num * potential(sys, qs[end, :]))

    scalar_num = dot(pseudosp.trial, num * pseudosp.trial)
    scalar_num_E = dot(pseudosp.trial, num_E * pseudosp.trial)
    scalar_denom = dot(sp.trial, denom * sp.trial)

    [scalar_num/scalar_denom,
     scalar_num_E/scalar_denom]
end

"""
Monte Carlo solution for an arbitrary PIGS system.
"""
struct PigsSampling <: Sampling
    "Pseudo-partition function."
    Z::Float64
    "Pseudo-partition function standard error."
    Z_err::Float64
    "Energy."
    E::Float64
    "Energy standard error."
    E_err::Float64
end

"""
    PigsSampling(sys::System, beta::Float64, P::Int, num_samples::Int; sampling_sys::Maybe{System}=nothing, progress_output::IO=stderr)

Calculate the solution for `sys` at `beta` with `P` links, `num_samples` random
samples, and a uniform (in space and surfaces) trial wavefunction.

If `sampling_sys` is provided, it is used for sampling. Otherwise, sampling
defaults to the simplified diagonal subsystem of `sys`.

The progress meter is written to `progress_output`.
"""
function PigsSampling(sys::System, beta::Float64, P::Int, num_samples::Int; sampling_sys::Maybe{System}=nothing, progress_output::IO=stderr)
    sys_diag = diag(sys)
    sys_diag_simple = simplify(sys_diag)
    pseudosp = PigsSamplingParameters(sys_diag_simple, beta, P)

    sampling_sys === nothing && (sampling_sys = sys_diag_simple)
    sp = PigsSamplingParameters(sampling_sys, beta, P)

    samples = zeros(Float64, 2, num_samples)
    meter = Progress(num_samples, output=progress_output)
    for n in ProgressWrapper(1:num_samples, meter)
        samples[:, n] .= get_sample_pigs(sys, pseudosp, sp)
    end

    simple = PigsAnalytical(sampling_sys, beta)
    normalization = simple.Z

    f_E(samples, samples_E) =
        samples_E ./ samples

    Z = mean(samples[1, :]) * normalization
    Z_err = std(samples[1, :]) / sqrt(num_samples) * normalization
    E, E_err = jackknife(f_E, [samples[n, :] for n in 1:size(samples, 1)]...)

    PigsSampling(Z, Z_err, E, E_err)
end
