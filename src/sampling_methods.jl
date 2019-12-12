# Methods for sampling from mixture distributions.

function loop_samples(f::Function, sp::AbstractSamplingParameters{S_,M,P}, num_samples::Int, progress_output::IO) where {S_,M,P}
    qs = zeros(num_beads(sp), M)
    meter = Progress(num_samples, output=progress_output)
    for n in ProgressWrapper(1:num_samples, meter)
        s_ = sample(1:S_, sp.weights)
        for m in 1:M
            qs[:, m] .= rand(sp.mvns[m, s_])
        end
        f(n, qs)
    end
end
