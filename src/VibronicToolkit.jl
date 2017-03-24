module VibronicToolkit

using Distributions: MvNormal, sample
import JSON
using ProgressMeter: @showprogress
using StatsBase: WeightVec

export
    Beta,
    System,

    Simple,

    Analytical,
    SumOverStates,
    Trotter,
    Sampling

"""
Inverse temperature.
"""
typealias Beta Float64
Beta(x::Float64) = x > 0.0 ? x : throw(DomainError())

"""
Solution for vibronic problem.
"""
abstract Solution

include("system.jl")
include("operators.jl")

include("simple.jl")
include("analytical.jl")
include("sos.jl")
include("trotter.jl")
include("sampling.jl")

end
