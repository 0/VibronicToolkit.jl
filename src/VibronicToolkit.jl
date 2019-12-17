module VibronicToolkit

import LinearAlgebra: diag, isdiag
using LinearAlgebra: diagm, Diagonal, dot, eigen, eigvals, I, issymmetric, svd,
                     Symmetric, tr
using Printf: @printf
using Statistics: mean, std, var

using Distributions: MvNormal, sample
import JSON
using ProgressMeter: Progress, ProgressWrapper
using Qutilities: ptrace
using Sobol: next!, SobolSeq
using SpecialFunctions: erfcinv
using StatsBase: Weights

# We re-export diag and isdiag from LinearAlgebra.
export
    Splitting,
    h0_V,
    p2_U,

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

    make_sampling_method,

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
Trotter splitting variant.
"""
@enum Splitting begin
    h0_V
    p2_U
end
@doc "Harmonic oscillators with low-order diagonal terms." h0_V
@doc "Kinetic energy." p2_U

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
include("spsa.jl")

include("analytical.jl")
include("sos.jl")
include("trotter.jl")

include("sampling.jl")

include("pes.jl")
include("density.jl")

include("iterative_decomposition.jl")

end
