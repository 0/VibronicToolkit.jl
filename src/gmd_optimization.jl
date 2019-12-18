# Optimization of GMD parameters using Monte Carlo at finite temperature.

"""
    get_sample_opt(sys::System{S,M}, pseudosp::SamplingParameters{S,M,P}, sp::SamplingParameters{S_,M,P}, qs::Matrix{Float64})

Compute a sample for `sys` using `pseudosp`, `sp`, and the points `qs`.
"""
function get_sample_opt(sys::System{S,M}, pseudosp::SamplingParameters{S,M,P}, sp::SamplingParameters{S_,M,P}, qs::Matrix{Float64}) where {S,S_,M,P}
    # Calculate the numerator and denominator matrices.
    num = Matrix{Float64}(I, S, S)
    denom = Matrix{Float64}(I, S_, S_)

    for i1 in 1:P
        # Next index in cyclic path.
        i2 = mod(i1+1, 1:P)

        # Numerator.
        FM_num, scaling = sampling_matrix_free_particle(pseudosp, qs[i1, :], qs[i2, :])
        _, IM = sampling_matrix_interaction(sys, pseudosp, qs[i1, :])
        num *= IM * FM_num

        # Denominator.
        FM_denom, _ = sampling_matrix_free_particle(sp, qs[i1, :], qs[i2, :]; scaling=scaling)
        denom *= FM_denom
    end

    trace_num = tr(num)
    trace_denom = tr(denom)

    [trace_num/trace_denom,
     abs(trace_num/trace_denom),
     log(abs(trace_num/trace_denom))]
end

"""
Sampling result for GMD optimization.
"""
struct SamplingGmdOptimization
    "Relative entropy."
    relent::Float64
    "Relative entropy standard error."
    relent_err::Float64
end

"""
    SamplingGmdOptimization(sys::System, beta::Float64, P::Int, sm::SamplingMethod; sampling_sys::System)

Calculate optimization quantities for `sys` at `beta` using `P` links, sampling
method `sm`, and sampling system `sampling_sys`.
"""
function SamplingGmdOptimization(sys::System, beta::Float64, P::Int, sm::SamplingMethod; sampling_sys::System)
    sys_diag = diag(sys)
    sys_diag_simple = simplify(sys_diag)
    pseudosp = SamplingParameters(sys_diag_simple, beta, P)
    sp = SamplingParameters(sampling_sys, beta, P)

    result = loop_samples(sp, sm, [:Z, :Z_abs, :Z_abs_log], devnull) do qs
        get_sample_opt(sys, pseudosp, sp, qs)
    end

    relent, relent_err = analyze(result, :Z_abs, :Z_abs_log) do sample_Z_abs, sample_Z_abs_log
        # The normalization cancels, so we don't worry about it.
        log(sample_Z_abs) - sample_Z_abs_log
    end

    SamplingGmdOptimization(relent, relent_err)
end

"""
Sampling system parameters used as inputs to the loss function.
"""
struct SpsaInputGmdOptimization{S,M} <: SpsaInput
    "Energy offsets (S, S)."
    energy::Matrix{Float64}
    "Linear prefactors (M, S, S)."
    lin::Array{Float64,3}
end

"""
    SpsaInputGmdOptimization(S::Int, M::Int)

Create zero-filled parameters for `S` surfaces and `M` modes.
"""
function SpsaInputGmdOptimization(S::Int, M::Int)
    energy = zeros(Float64, S, S)
    lin = zeros(Float64, M, S, S)
    SpsaInputGmdOptimization{S,M}(energy, lin)
end

"""
    SpsaInputGmdOptimization(sys::System)

Extract parameters from `sys`.
"""
function SpsaInputGmdOptimization(sys::System{S,M}) where {S,M}
    SpsaInputGmdOptimization{S,M}(sys.coef[0], sys.coef[1])
end

function Base.:+(si1::SpsaInputGmdOptimization{S,M}, si2::SpsaInputGmdOptimization{S,M}) where {S,M}
    SpsaInputGmdOptimization{S,M}(si1.energy + si2.energy, si1.lin + si2.lin)
end

function Base.:*(x::Float64, si::SpsaInputGmdOptimization{S,M}) where {S,M}
    SpsaInputGmdOptimization{S,M}(x * si.energy, x * si.lin)
end

function Base.rand(si::SpsaInputGmdOptimization{S,M}) where {S,M}
    energy = zero(si.energy)
    lin = zero(si.lin)
    for s in 1:S
        energy[s, s] = bernoulli()
        for m in 1:M
            lin[m, s, s] = bernoulli()
        end
    end
    SpsaInputGmdOptimization{S,M}(energy, lin)
end

function DiagonalSystem(freq::Matrix{Float64}, si::SpsaInputGmdOptimization)
    DiagonalSystem(freq, [si.energy, si.lin])
end

"""
Parameters for GMD optimization.
"""
struct GmdOptimizationParameters{P}
    "Reciprocal temperature."
    beta::Float64
    "Sampling method."
    sm::SamplingMethod

    "Lower bound on desired loss function values."
    loss_bound_lower::Float64
    "Upper bound on desired loss function values."
    loss_bound_upper::Float64

    "Surface repulsion stiffness."
    repulsion_stiffness::Float64
    "Surface repulsion penalty coefficient."
    repulsion_coefficient::Float64
    "Weight imbalance penalty coefficient."
    weight_imbalance_coefficient::Float64
end

"""
    loss(gop::GmdOptimizationParameters, sys::System{S,M}, sampling_sys::DiagonalSystem{S_,M})

Compute the loss function for `sys` using `sampling_sys` and `gop`.
"""
function loss(gop::GmdOptimizationParameters{P}, sys::System{S,M}, sampling_sys::DiagonalSystem{S_,M}) where {P,S,M,S_}
    sampling = SamplingGmdOptimization(sys, gop.beta, P, gop.sm; sampling_sys=sampling_sys)
    l = sampling.relent
    err = sampling.relent_err

    for s_1 in 1:(S_-1)
        for s_2 in (s_1+1):S_
            dist2 = sum((sampling_sys.ds[:, s_1] .- sampling_sys.ds[:, s_2]).^2)
            weight = exp(-0.5 * dist2 / gop.repulsion_stiffness^2)
            l += gop.repulsion_coefficient * 2 * weight / (S_ * (S_ - 1))
        end
    end

    weight_imbalance = (1.0 - shannon_norm(weights(sampling_sys, gop.beta)))
    l += gop.weight_imbalance_coefficient * weight_imbalance

    l, err
end

"""
Optimized sampling GMD results.
"""
abstract type AbstractGmdOptimization end

DiagonalSystem(opt::AbstractGmdOptimization) = DiagonalSystem(opt.freq, opt.si)

"""
Optimized sampling GMD results without deformation.
"""
struct GmdOptimization <: AbstractGmdOptimization
    "Frequencies (M, S)."
    freq::Matrix{Float64}
    "System parameters."
    si::SpsaInputGmdOptimization
end

"""
    GmdOptimization(sys::System{S,M}, gop::GmdOptimizationParameters, num_iter::Int, spsa_a::Float64, start_sys::DiagonalSystem{S_,M}; progress_output::IO=stderr, data_output::IO=devnull)

Optimize `start_sys` to better match `sys` using `gop`, `num_iter`, and
`spsa_a`.

The progress meter is written to `progress_output`. Loss function data is
written to `data_output`.
"""
function GmdOptimization(sys::System{S,M}, gop::GmdOptimizationParameters, num_iter::Int, spsa_a::Float64, start_sys::DiagonalSystem{S_,M}; progress_output::IO=stderr, data_output::IO=devnull) where {S,M,S_}
    si = SpsaInputGmdOptimization(start_sys)
    freq = start_sys.freq

    println(data_output, "# iter loss err")
    flush(data_output)

    function loss_(si_::SpsaInputGmdOptimization)
        loss(gop, sys, DiagonalSystem(freq, si_))
    end

    sp, _, _ = SpsaParameters(loss_, num_iter, spsa_a, si)

    meter = Progress(num_iter, output=progress_output)
    si = spsa(sp, si) do k, si_
        l, err = loss_(si_)
        println(data_output, "$(k) $(l) $(err)")
        flush(data_output)
        update!(meter, k)
    end
    finish!(meter)

    GmdOptimization(freq, si)
end

"""
    step_nu(gop::GmdOptimizationParameters, sys::System, sampling_sys::System, nu_prev::Float64, dnu_prev::Float64)

Increase `nu_prev` by an amount that satisfies the bounds in `gop` for `sys`
and `sampling_sys`, using `dnu_prev` as a guideline.
"""
function step_nu(gop::GmdOptimizationParameters, sys::System, sampling_sys::System, nu_prev::Float64, dnu_prev::Float64)
    dnu_min = 1e-3
    if nu_prev > 1.0 - 1e-1 - 1e-2
        dnu_max = 1e-2
    else
        dnu_max = 1e-1
    end
    dnu = min(dnu_max, 2dnu_prev)
    for i in 1:8
        nu = min(1.0, nu_prev + dnu)
        l, _ = loss(gop, nu * sys, sampling_sys)
        if l < gop.loss_bound_lower
            dnu_min = dnu
        elseif l > gop.loss_bound_upper
            dnu_max = dnu
        else
            break
        end
        dnu = 0.5 * (dnu_min + dnu_max)
    end
    nu = min(1.0, nu_prev + dnu)
    nu, nu - nu_prev
end

"""
Optimized sampling GMD results with deformation.
"""
struct GmdOptimizationDeformation <: AbstractGmdOptimization
    "Frequencies (M, S)."
    freq::Matrix{Float64}
    "System parameters."
    si::SpsaInputGmdOptimization
end

"""
    GmdOptimizationDeformation(sys::System, gop::GmdOptimizationParameters, num_iter::Int, spsa_a::Float64, S_::Int; progress_output::IO=stderr, data_output::IO=devnull)

Optimize an empty sampling system with `S_` surfaces to match `sys` using
`gop`, `num_iter`, and `spsa_a`.

The progress meter is written to `progress_output`. Loss function data is
written to `data_output`.
"""
function GmdOptimizationDeformation(sys::System{S,M}, gop::GmdOptimizationParameters, num_iter::Int, spsa_a::Float64, S_::Int; progress_output::IO=stderr, data_output::IO=devnull) where {S,M}
    si = SpsaInputGmdOptimization(S_, M)
    freq = zeros(Float64, M, S_)
    for s_ in 1:S_
        freq[:, s_] .= sys.freq[:, 1]
    end

    num_walkers = Threads.nthreads()

    println(data_output, "# step nu iter loss err")
    flush(data_output)

    step = 0
    nu = 0.0
    dnu = 1e-1

    meter = Progress(10_000, output=progress_output)
    while nu < 1.0
        step += 1
        nu, dnu = step_nu(gop, sys, DiagonalSystem(freq, si), nu, dnu)
        sys_nu = nu * sys

        function loss_(si_::SpsaInputGmdOptimization)
            loss(gop, sys_nu, DiagonalSystem(freq, si_))
        end

        losses = fill(Inf, num_iter + 1, num_walkers)
        errs = fill(Inf, num_iter + 1, num_walkers)
        sis = Array{SpsaInputGmdOptimization}(undef, num_iter + 1, num_walkers)

        sp, l, err = SpsaParameters(loss_, num_iter, spsa_a, si)

        losses[1, :] .= l
        errs[1, :] .= err
        sis[1, :] .= Ref(si)

        Threads.@threads for _ in 1:num_walkers
            thridx = Threads.threadid()

            spsa(sp, si) do k, si_
                losses[k+1, thridx], errs[k+1, thridx] = loss_(si_)
                sis[k+1, thridx] = si_
            end
        end

        I = argmin([l + err for (l, err) in zip(losses, errs)])
        best_kp1 = I[1]
        best_walker = I[2]
        si = sis[I]

        for kp1 in 1:best_kp1
            l = losses[kp1, best_walker]
            err = errs[kp1, best_walker]
            println(data_output, "$(step) $(nu) $(kp1-1) $(l) $(err)")
            flush(data_output)
        end

        update!(meter, Int(floor(nu * meter.n)))
    end
    finish!(meter)

    GmdOptimizationDeformation(freq, si)
end
