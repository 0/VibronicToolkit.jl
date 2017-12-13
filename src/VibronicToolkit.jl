module VibronicToolkit

using Distributions: MvNormal, sample
import JSON
using ProgressMeter: @showprogress
using StatsBase: Weights

export
    System,

    simplify,

    Analytical,
    SumOverStates,
    Trotter,
    SamplingFiniteDifference,
    SamplingPrimitiveThermodynamic

"""
Solution for vibronic problem.
"""
abstract type Solution end

include("system.jl")
include("operators.jl")

include("jackknife.jl")

include("analytical.jl")
include("sos.jl")
include("trotter.jl")

include("sampling.jl")
include("sampling_finite_difference.jl")
include("sampling_primitive_thermodynamic.jl")

end
