#!/usr/bin/env julia

using VibronicToolkit

using ArgParse

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
    "--pigs"
        help = "use PIGS instead of finite temperature"
        action = :store_true
    "--conf"
        metavar = "FILE"
        help = "path to config file"
        required = true
    "--beta"
        metavar = "T"
        help = "reciprocal temperature"
        arg_type = Float64
        required = true
    "--dbeta"
        metavar = "T"
        help = "finite difference step (in units of beta)"
        arg_type = Float64
    "--num-links"
        metavar = "P"
        help = "number of Trotter links"
        arg_type = Int
        required = true
    "--num-samples"
        metavar = "N"
        help = "number of samples"
        arg_type = Int
        required = true
    "--sampling-conf"
        metavar = "FILE"
        help = "path to sampling config file"
    "--quiet"
        help = "suppress progress meter"
        action = :store_true
end
c = parse_args(ARGS, s, as_symbols=true)

sys = read(c[:conf], DenseSystem)
beta = c[:beta]
if c[:dbeta] !== nothing
    dbeta = c[:dbeta] * beta
else
    dbeta = nothing
end
P = c[:num_links]
num_samples = c[:num_samples]
if c[:sampling_conf] !== nothing
    sampling_sys = read(c[:sampling_conf], DiagonalSystem)
else
    sampling_sys = nothing
end
if c[:quiet]
    progress_output=devnull
else
    progress_output=stderr
end

if c[:pigs]
    trial = UniformTrialWavefunction(sys)

    if dbeta !== nothing
        error("Finite difference not supported for PIGS")
    else
        sampling = PigsSampling(sys, trial, beta, P, num_samples; sampling_sys=sampling_sys, progress_output=progress_output)
    end

    println("Z: $(sampling.Z) ± $(sampling.Z_err)")
    println("E: $(sampling.E) ± $(sampling.E_err)")
    println("SvN: $(sampling.SvN) ± $(sampling.SvN_err)")
    println("S2: $(sampling.S2) ± $(sampling.S2_err)")
else
    if dbeta !== nothing
        sampling = SamplingFiniteDifference(sys, beta, dbeta, P, num_samples; sampling_sys=sampling_sys, progress_output=progress_output)
    else
        sampling = SamplingPrimitiveThermodynamic(sys, beta, P, num_samples; sampling_sys=sampling_sys, progress_output=progress_output)
    end

    println("Z: $(sampling.Z) ± $(sampling.Z_err)")
    println("E: $(sampling.E) ± $(sampling.E_err)")
    println("Cv: $(sampling.Cv) ± $(sampling.Cv_err)")
end
