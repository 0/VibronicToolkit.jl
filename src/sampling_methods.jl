# Methods for sampling from mixture distributions.

"""
Method of obtaining samples from surfaces.
"""
abstract type SamplingMethod end

"""
Result obtained by sampling from surfaces.
"""
abstract type SamplingResult end

function Base.getindex(sr::SamplingResult, mean_label::Symbol)
    idx = findfirst(isequal(mean_label), sr.mean_labels)
    sr.samples[idx, :]
end

"""
Choose a random surface, then draw a random sample from it.
"""
struct SamplingMethodRandom <: SamplingMethod
    "Number of samples drawn from all surfaces."
    num_samples::Int

    function SamplingMethodRandom(num_samples::Int)
        num_samples >= 1 || throw(DomainError(num_samples, "At least one sample."))
        new(num_samples)
    end
end

"""
Result for `SamplingMethodRandom`.
"""
struct SamplingResultRandom <: SamplingResult
    mean_labels::Vector{Symbol}
    samples::Array{Any,2}
end

function loop_samples(f::Function, sp::AbstractSamplingParameters{S_,M,P}, sm::SamplingMethodRandom, mean_labels::Vector{Symbol}, progress_output::IO) where {S_,M,P}
    qs = zeros(num_beads(sp), M)
    samples = Array{Any}(undef, length(mean_labels), sm.num_samples)

    meter = Progress(sm.num_samples, output=progress_output)
    for n in ProgressWrapper(1:sm.num_samples, meter)
        s_ = sample(1:S_, sp.weights)
        for m in 1:M
            qs[:, m] .= rand(sp.mvns[m, s_])
        end
        samples[:, n] .= f(qs)
    end

    SamplingResultRandom(mean_labels, samples)
end

function analyze(f::Function, sr::SamplingResultRandom, mean_labels::Symbol...)
    data = [sr[mean_label] for mean_label in mean_labels]
    jackknife(f, data...)
end

"""
Choose the best surface, then draw a random sample from it.
"""
struct SamplingMethodDeterministic <: SamplingMethod
    "Number of samples drawn from all surfaces."
    num_samples::Int
    "Number of bootstrap samples for each surface."
    num_boot::Int

    function SamplingMethodDeterministic(num_samples::Int, num_boot::Int)
        num_samples >= 1 || throw(DomainError(num_samples, "At least one sample."))
        num_boot >= 2 || throw(DomainError(num_boot, "At least two bootstrap samples."))
        new(num_samples, num_boot)
    end
end

"""
Result for `SamplingMethodDeterministic`.
"""
struct SamplingResultDeterministic <: SamplingResult
    mean_labels::Vector{Symbol}
    samples::Array{Any,2}
    weights::Weights
    surface_labels::Vector{Int}
end

function loop_samples(f::Function, sp::AbstractSamplingParameters{S_,M,P}, sm::SamplingMethodDeterministic, mean_labels::Vector{Symbol}, progress_output::IO) where {S_,M,P}
    sm.num_samples >= S_ * sm.num_boot || throw(DomainError(sm, "Not enough samples."))
    num_samples_boot = S_ * sm.num_boot

    !isempty(mean_labels) && mean_labels[1] == :Z || throw(DomainError(mean_labels, "Z must be first."))

    qs = zeros(num_beads(sp), M)
    samples = Array{Any}(undef, length(mean_labels), sm.num_samples)
    surface_labels = zeros(Int, sm.num_samples)

    Ns = zeros(Int, S_)
    xs_boot = zeros(Float64, sm.num_boot, S_)
    means_boot = zeros(Float64, S_)
    Ts = zeros(Float64, S_)
    Us = zeros(Float64, S_)
    err2s = zeros(Float64, S_)
    deltas = zeros(Float64, S_)

    meter = Progress(sm.num_samples, output=progress_output)
    for n in ProgressWrapper(1:sm.num_samples, meter)
        if n > num_samples_boot
            s_ = argmax(deltas)
        else
            # Bootstrap.
            s_ = 1 + div(n - 1, sm.num_boot)
        end

        Ns[s_] += 1
        surface_labels[n] = s_
        for m in 1:M
            qs[:, m] .= rand(sp.mvns[m, s_])
        end
        samples[:, n] .= f(qs)
        x = samples[1, n]

        if n > num_samples_boot
            delta_x = x - means_boot[s_]
            Ts[s_] += delta_x
            Us[s_] += (Ns[s_] * delta_x - Ts[s_])^2 / (Ns[s_] * (Ns[s_] + 1))
        else
            # Bootstrap.
            xs_boot[Ns[s_], s_] = x

            if Ns[s_] == sm.num_boot
                # End of bootstrap for current surface.
                means_boot[s_] = mean(xs_boot[:, s_])
                Us[s_] = sum((xs_boot[:, s_] .- means_boot[s_]).^2)
            end
        end

        err2s[s_] = Us[s_] / (Ns[s_] * (Ns[s_] + 1))
        deltas[s_] = sp.weights[s_]^2 * err2s[s_] / (Ns[s_] + 1)
    end

    SamplingResultDeterministic(mean_labels, samples, sp.weights, surface_labels)
end

function analyze(f::Function, sr::SamplingResultDeterministic, mean_labels::Symbol...)
    data = [sr[mean_label] for mean_label in mean_labels]
    jackknife(f, sr.weights, sr.surface_labels, data...)
end

"""
Choose the best surface, then draw a quasi-random sample from it.
"""
struct SamplingMethodQuasirandom <: SamplingMethod
    "Number of samples drawn from all surfaces."
    num_samples::Int
    "Number of bootstrap samples for each surface."
    num_boot::Int
    "Number of randomized sequences."
    num_seqs::Int

    function SamplingMethodQuasirandom(num_samples::Int, num_boot::Int, num_seqs::Int)
        num_samples >= 1 || throw(DomainError(num_samples, "At least one sample."))
        num_boot >= 1 || throw(DomainError(num_boot, "At least one bootstrap sample."))
        num_seqs >= 2 || throw(DomainError(num_seqs, "At least two sequences."))
        new(num_samples, num_boot, num_seqs)
    end
end

"""
Result for `SamplingMethodQuasirandom`.
"""
struct SamplingResultQuasirandom <: SamplingResult
    mean_labels::Vector{Symbol}
    samples::Array{Any,2}
    weights::Weights
    surface_labels::Vector{Int}
end

function loop_samples(f::Function, sp::AbstractSamplingParameters{S_,M,P}, sm::SamplingMethodQuasirandom, mean_labels::Vector{Symbol}, progress_output::IO) where {S_,M,P}
    sm.num_samples >= S_ * sm.num_boot || throw(DomainError(sm, "Not enough samples."))
    num_samples_boot = S_ * sm.num_boot

    !isempty(mean_labels) && mean_labels[1] == :Z || throw(DomainError(mean_labels, "Z must be first."))

    Q = num_beads(sp)

    seqs = [SobolSeq(Q*M) for s in 1:S_]
    shifts = rand(Q, M, sm.num_seqs, S_)

    eigvals = Array{Float64}(undef, Q, M, S_)
    eigvecs = Array{Float64}(undef, Q, Q, M, S_)
    for s_ in 1:S_
        for m in 1:M
            F = eigen(sp.precs[m, s_])
            eigvals[:, m, s_] .= F.values
            eigvecs[:, :, m, s_] .= F.vectors
        end
    end

    Ns = zeros(Int, S_)
    Ts = Array{Any}(undef, length(mean_labels), sm.num_seqs, S_)
    pre_means = Array{Any}(undef, length(mean_labels), sm.num_seqs, S_)
    means = zeros(Float64, S_)
    err2s = zeros(Float64, S_)
    deltas = zeros(Float64, S_)

    u = zeros(Q, M)
    ut = zeros(Q, M)
    qs = zeros(Q, M)

    meter = Progress(sm.num_samples, output=progress_output)
    for n in ProgressWrapper(1:sm.num_samples, meter)
        if n > num_samples_boot
            s_ = argmax(deltas)
        else
            # Bootstrap.
            s_ = 1 + div(n - 1, sm.num_boot)
        end

        Ns[s_] += 1

        u .= reshape(next!(seqs[s_]), Q, M)
        for t in 1:sm.num_seqs
            for m in 1:M
                for j in 1:Q
                    phi = (u[j, m] + shifts[j, m, t, s_]) % 1
                    ut[j, m] = gauss_cdf_inv(phi) / sqrt(eigvals[j, m, s_])
                end
                qs[:, m] .= eigvecs[:, :, m, s_] * ut[:, m] + sp.means[m, s_]
            end

            if Ns[s_] == 1
                # Initialize.
                Ts[:, t, s_] .= f(qs)
            else
                Ts[:, t, s_] .+= f(qs)
            end
            pre_means[:, t, s_] .= Ts[:, t, s_] ./ Ns[s_]
        end

        means[s_] = mean(pre_means[1, :, s_])
        err2s[s_] = var(pre_means[1, :, s_]) / sm.num_seqs
        deltas[s_] = sp.weights[s_]^2 * err2s[s_] / (Ns[s_] + 1)
    end

    samples = reshape(pre_means, length(mean_labels), sm.num_seqs * S_)
    surface_labels = [s_ for s_ in 1:S_ for _ in 1:sm.num_seqs]

    SamplingResultQuasirandom(mean_labels, samples, sp.weights, surface_labels)
end

function analyze(f::Function, sr::SamplingResultQuasirandom, mean_labels::Symbol...)
    data = [sr[mean_label] for mean_label in mean_labels]
    jackknife(f, sr.weights, sr.surface_labels, data...)
end

function make_sampling_method(num_samples, num_boot, num_seqs)
    if !isnothing(num_seqs)
        SamplingMethodQuasirandom(num_samples, num_boot, num_seqs)
    elseif !isnothing(num_boot)
        SamplingMethodDeterministic(num_samples, num_boot)
    else
        SamplingMethodRandom(num_samples)
    end
end
