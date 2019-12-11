# Monte Carlo solution for finite temperature with primitive thermodynamic
# estimators.

"""
    get_sample_pt(sys::System{S,M}, pseudosp::SamplingParameters{S,M,P}, sp::SamplingParameters{S_,M,P}, qs::Matrix{Float64})

Compute a sample for `sys` using `pseudosp`, `sp`, and the points `qs`.
"""
function get_sample_pt(sys::System{S,M}, pseudosp::SamplingParameters{S,M,P}, sp::SamplingParameters{S_,M,P}, qs::Matrix{Float64}) where {S,S_,M,P}
    # Calculate the numerator, energy, heat capacity, and denominator matrices.
    num = Matrix{Float64}(I, S, S)
    num_Es = [Matrix{Float64}(I, S, S) for _ in 1:P]
    num_Cvs = [Matrix{Float64}(I, S, S) for _ in 1:P]
    denom = Matrix{Float64}(I, S_, S_)

    for i1 in 1:P
        # Next index in cyclic path.
        i2 = mod(i1+1, 1:P)

        # Numerator.
        FM_num, scaling = sampling_matrix_free_particle(pseudosp, qs[i1, :], qs[i2, :])
        EM1, IM = sampling_matrix_interaction(sys, pseudosp, qs[i1, :])
        EM2 = sampling_matrix_energy(pseudosp, qs[i1, :], qs[i2, :])
        EM = EM1 - EM2
        if i1 == 1
            CM1 = sampling_matrix_heat_capacity(pseudosp, qs[i1, :], qs[i2, :])
            CM = EM1 * EM + CM1 - EM * EM2
        end
        IFM = IM * FM_num
        IEFM = IM * EM * FM_num
        num *= IM * FM_num
        for i in 1:P
            if i == i1
                num_Es[i] *= IEFM
                if i1 == 1
                    num_Cvs[i] *= IM * CM * FM_num
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

        # Denominator.
        FM_denom, _ = sampling_matrix_free_particle(sp, qs[i1, :], qs[i2, :]; scaling=scaling)
        denom *= FM_denom
    end

    trace_num = tr(num)
    trace_num_E = mean(tr(num_E) for num_E in num_Es)
    trace_num_Cv = mean(tr(num_Cv) for num_Cv in num_Cvs)
    trace_denom = tr(denom)

    [trace_num/trace_denom,
     trace_num_E/trace_denom,
     trace_num_Cv/trace_denom]
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
    SamplingPrimitiveThermodynamic(sys::System, beta::Float64, P::Int, num_samples::Int; sampling_sys::Maybe{System}=nothing, progress_output::IO=stderr)

Calculate the solution for `sys` at `beta` using `P` links and `num_samples`
random samples.

If `sampling_sys` is provided, it is used for sampling. Otherwise, sampling
defaults to the simplified diagonal subsystem of `sys`.

The progress meter is written to `progress_output`.
"""
function SamplingPrimitiveThermodynamic(sys::System, beta::Float64, P::Int, num_samples::Int; sampling_sys::Maybe{System}=nothing, progress_output::IO=stderr)
    sys_diag = diag(sys)
    sys_diag_simple = simplify(sys_diag)
    pseudosp = SamplingParameters(sys_diag_simple, beta, P)

    isnothing(sampling_sys) && (sampling_sys = sys_diag_simple)
    sp = SamplingParameters(sampling_sys, beta, P)

    samples_Z = zeros(Float64, num_samples)
    samples_E = zeros(Float64, num_samples)
    samples_Cv = zeros(Float64, num_samples)

    loop_samples(sp, num_samples, progress_output) do n, qs
        sample = get_sample_pt(sys, pseudosp, sp, qs)
        samples_Z[n] = sample[1]
        samples_E[n] = sample[2]
        samples_Cv[n] = sample[3]
    end

    simple = Analytical(sampling_sys, beta)
    normalization = simple.Z

    Z, Z_err = jackknife(samples_Z) do sample_Z
        sample_Z * normalization
    end

    E, E_err = jackknife(samples_Z, samples_E) do sample_Z, sample_E
        sample_E / sample_Z
    end

    Cv, Cv_err = jackknife(samples_Z, samples_E, samples_Cv) do sample_Z, sample_E, sample_Cv
        (sample_Cv / sample_Z - (sample_E / sample_Z)^2) * beta^2
    end

    SamplingPrimitiveThermodynamic(Z, Z_err, E, E_err, Cv, Cv_err)
end
