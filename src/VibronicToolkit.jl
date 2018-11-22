module VibronicToolkit

import LinearAlgebra: diag, isdiag
using LinearAlgebra: diagm, Diagonal, dot, eigen, eigvals, I, issymmetric, svd,
                     Symmetric, tr
using Printf: @printf
using Statistics: mean, std

using Distributions: MvNormal, sample
import JSON
using ProgressMeter: Progress, ProgressWrapper
using Qutilities: ptrace
using StatsBase: Weights

# We re-export diag and isdiag from LinearAlgebra.
export
    DegeneracyException,
    SurfaceCouplingException,

    System,
    DenseSystem,
    DiagonalSystem,

    diag,
    isdiag,
    simplify,
    weights,

    TrialWavefunction,
    UniformTrialWavefunction,

    Analytical,
    SumOverStates,
    Trotter,
    SamplingFiniteDifference,
    SamplingPrimitiveThermodynamic,

    PigsAnalytical,
    PigsSumOverStates,
    PigsTrotter,
    PigsSampling,

    GroundPes,
    path_mean_std,
    DiagonalDensity,
    PigsDiagonalDensity,

    IterativeDecomposition

"""
Solution for vibronic problem.
"""
abstract type Solution end

include("utilities.jl")
include("entanglement.jl")

include("system.jl")
include("operators.jl")
include("trial.jl")

include("jackknife.jl")

include("analytical.jl")
include("sos.jl")
include("trotter.jl")

include("sampling.jl")
include("sampling_finite_difference.jl")
include("sampling_pigs.jl")
include("sampling_primitive_thermodynamic.jl")

include("pes.jl")
include("density.jl")

include("iterative_decomposition.jl")

end
