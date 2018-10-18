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
    num_L = Matrix{Float64}(I, S, S)
    num_R = Matrix{Float64}(I, S, S)
    denom = Matrix{Float64}(I, S_, S_)

    middle_bead = P/2 + 1

    for i1 in 1:P
        # Next index in open path.
        i2 = i1+1

        # Numerator.
        FM_num, scaling = sampling_matrix_free_particle(pseudosp, qs[i1, :], qs[i2, :])
        _, IM_L = sampling_matrix_interaction(sys, pseudosp, qs[i1, :])
        _, IM_R = sampling_matrix_interaction(sys, pseudosp, qs[i2, :])

        if i1 < middle_bead
            num_L *= sqrt(IM_L) * FM_num * sqrt(IM_R)
        else
            num_R *= sqrt(IM_L) * FM_num * sqrt(IM_R)
        end

        # Denominator.
        FM_denom, _ = sampling_matrix_free_particle(sp, qs[i1, :], qs[i2, :]; scaling=scaling)
        denom *= FM_denom
    end

    num = num_L * num_R

    # Average the energy at the two ends.
    num_E = 0.5 * (potential(sys, qs[1, :]) * num + num * potential(sys, qs[end, :]))

    trial_vec_num_L = trial_spatial(pseudosp.trial, qs[1, :])
    trial_vec_num_R = trial_spatial(pseudosp.trial, qs[end, :])
    trial_vec_denom_L = trial_spatial(sp.trial, qs[1, :])
    trial_vec_denom_R = trial_spatial(sp.trial, qs[end, :])

    scalar_num = dot(trial_vec_num_L, num * trial_vec_num_R)
    scalar_num_E = dot(trial_vec_num_L, num_E * trial_vec_num_R)
    scalar_denom = dot(trial_vec_denom_L, denom * trial_vec_denom_R)

    reduced_density_matrix = (num_R * trial_vec_num_R) * (trial_vec_num_L' * num_L)
    # Symmetrize.
    reduced_density_matrix = 0.5 * (reduced_density_matrix + reduced_density_matrix')

    [scalar_num/scalar_denom,
     scalar_num_E/scalar_denom,
     reduced_density_matrix/scalar_denom]
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
    "Von Neumann entanglement entropy."
    SvN::Float64
    "Von Neumann entanglement entropy standard error."
    SvN_err::Float64
    "Order-2 Rényi entanglement entropy."
    S2::Float64
    "Order-2 Rényi entanglement entropy standard error."
    S2_err::Float64
end

"""
    PigsSampling(sys::System, trial::TrialWavefunction, beta::Float64, P::Int, num_samples::Int; sampling_sys::Maybe{System}=nothing, progress_output::IO=stderr)

Calculate the solution for `sys` with `trial` propagated by `beta` using `P`
links and `num_samples` random samples.

If `sampling_sys` is provided, it is used for sampling. Otherwise, sampling
defaults to the simplified diagonal subsystem of `sys`.

The progress meter is written to `progress_output`.
"""
function PigsSampling(sys::System, trial::TrialWavefunction, beta::Float64, P::Int, num_samples::Int; sampling_sys::Maybe{System}=nothing, sampling_trial::Maybe{TrialWavefunction}=nothing, progress_output::IO=stderr)
    sys_diag = diag(sys)
    sys_diag_simple = simplify(sys_diag)
    pseudosp = PigsSamplingParameters(sys_diag_simple, trial, beta, P)

    sampling_sys === nothing && (sampling_sys = sys_diag_simple)
    sampling_trial === nothing && (sampling_trial = trial)
    sp = PigsSamplingParameters(sampling_sys, sampling_trial, beta, P)

    samples = zeros(Float64, 2, num_samples)
    rhos = Matrix{Float64}[]
    meter = Progress(num_samples, output=progress_output)
    for n in ProgressWrapper(1:num_samples, meter)
        sample = get_sample_pigs(sys, pseudosp, sp)
        samples[:, n] .= sample[1:2]
        push!(rhos, sample[3])
    end

    simple = PigsAnalytical(sampling_sys, sampling_trial, beta)
    normalization = simple.Z

    f_E(sample, sample_E) =
        sample_E / sample
    f_SvN(sample, rho) =
        S_vn(rho / sample)
    f_S2(sample, rho) =
        S_renyi(rho / sample)

    Z = mean(samples[1, :]) * normalization
    Z_err = std(samples[1, :]) / sqrt(num_samples) * normalization
    E, E_err = jackknife(f_E, [samples[n, :] for n in 1:size(samples, 1)]...)
    SvN, SvN_err = jackknife(f_SvN, samples[1, :], rhos)
    S2, S2_err = jackknife(f_S2, samples[1, :], rhos)

    PigsSampling(Z, Z_err, E, E_err, SvN, SvN_err, S2, S2_err)
end
