module VibronicToolkit

using Distributions: MvNormal, sample
import JSON
using ProgressMeter: @showprogress
using StatsBase: WeightVec

export
    System,

    simplify,

    Analytical,
    SumOverStates,
    Trotter,
    Sampling

"""
Solution for vibronic problem.
"""
abstract Solution

include("system.jl")
include("operators.jl")

include("analytical.jl")
include("sos.jl")
include("trotter.jl")
include("sampling.jl")

end
