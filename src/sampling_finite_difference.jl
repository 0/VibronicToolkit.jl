# Monte Carlo solution with finite difference estimators.

"""
    get_sample_fd(sys::System{S,M}, sps::SamplingParameters{S,M,P}...)

Compute a sample for `sys` using `sps`.
"""
function get_sample_fd(sys::System{S,M}, sps::SamplingParameters{S,M,P}...) where {S,M,P}
    length(sps) >= 1 || throw(DomainError(length(sps), "At least one set of parameters."))

    # Choose a surface.
    s_ = sample(1:S, sps[1].weights)

    # Sample coordinates.
    qs = zeros(P, M)
    for m in 1:M
        qs[:, m] .= rand(sps[1].mvns[m, s_])
    end

    traces_num = Float64[]
    traces_denom = Float64[]

    scalings = Any[nothing for _ in 1:P]
    for sp in sps
        # Calculate the numerator and denominator matrices.
        num = Matrix{Float64}(I, S, S)
        denom = Matrix{Float64}(I, S, S)

        for i1 in 1:P
            # Next index in cyclic path.
            i2 = i1%P+1

            FM, scalings[i1] = sampling_matrix_free_particle(sp, qs[i1, :], qs[i2, :]; scaling=scalings[i1])
            _, IM = sampling_matrix_interaction(sys, sp, qs[i1, :])

            num *= IM * FM
            denom *= FM
        end

        push!(traces_num, tr(num))
        push!(traces_denom, tr(denom))
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
    SamplingFiniteDifference(sys::System, beta::Float64, dbeta::Float64, P::Int, num_samples::Int; progress_output::IO=stderr)

Calculate the solution for `sys` at `beta` with `P` links and `num_samples`
random samples, using finite difference step `dbeta`.

The progress meter is written to `progress_output`.
"""
function SamplingFiniteDifference(sys::System, beta::Float64, dbeta::Float64, P::Int, num_samples::Int; progress_output::IO=stderr)
    sys_diag = diag(sys)
    sys_diag_simple = simplify(sys_diag)
    sp = SamplingParameters(sys_diag_simple, beta, P)
    sp_m = SamplingParameters(sys_diag_simple, beta-dbeta, P)
    sp_p = SamplingParameters(sys_diag_simple, beta+dbeta, P)

    samples = zeros(Float64, 3, num_samples)
    meter = Progress(num_samples, output=progress_output)
    for n in ProgressWrapper(1:num_samples, meter)
        samples[:, n] .= get_sample_fd(sys, sp, sp_m, sp_p)
    end

    simple = Analytical(sys_diag_simple, beta)
    simple_m = Analytical(sys_diag_simple, beta-dbeta)
    simple_p = Analytical(sys_diag_simple, beta+dbeta)
    Zrat_m = simple.Z / simple_m.Z
    Zrat_p = simple.Z / simple_p.Z

    f_E(samples, samples_m, samples_p) =
        simple.E .+
        1.0/(2dbeta) * (Zrat_m * samples_m .- Zrat_p * samples_p) ./ samples
    f_Cv(samples, samples_m, samples_p) =
        simple.Cv .+
        (1.0/dbeta^2 * (Zrat_m * samples_m .+ Zrat_p * samples_p) ./ samples .-
         2.0/dbeta^2 .-
         (E - simple.E)^2) * beta^2

    Z = mean(samples[1, :])
    Z_err = std(samples[1, :]) / sqrt(num_samples)
    E, E_err = jackknife(f_E, [samples[n, :] for n in 1:size(samples, 1)]...)
    Cv, Cv_err = jackknife(f_Cv, [samples[n, :] for n in 1:size(samples, 1)]...)

    SamplingFiniteDifference(Z, Z_err, E, E_err, Cv, Cv_err)
end
