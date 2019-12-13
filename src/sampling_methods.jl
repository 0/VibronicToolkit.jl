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

function make_sampling_method(num_samples)
    SamplingMethodRandom(num_samples)
end
