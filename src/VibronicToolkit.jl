module VibronicToolkit

import LinearAlgebra: diag, isdiag
using LinearAlgebra: diagm, Diagonal, eigen, eigvals, I, issymmetric,
                     Symmetric, tr
using Printf: @printf
using Statistics: mean, std

using Distributions: MvNormal, sample
import JSON
using ProgressMeter: @showprogress
using StatsBase: Weights

# We re-export diag and isdiag from LinearAlgebra.
export
    SurfaceCouplingException,

    System,
    DenseSystem,
    DiagonalSystem,

    diag,
    isdiag,
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

include("utilities.jl")

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
