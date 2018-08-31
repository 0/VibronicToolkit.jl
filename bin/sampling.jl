#!/usr/bin/env julia

using VibronicToolkit

using ArgParse

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
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

if dbeta !== nothing
    sampling = SamplingFiniteDifference(sys, beta, dbeta, P, num_samples)
else
    sampling = SamplingPrimitiveThermodynamic(sys, beta, P, num_samples)
end

println("Z/Z: $(sampling.Z) ± $(sampling.Z_err)")
println("E: $(sampling.E) ± $(sampling.E_err)")
println("Cv: $(sampling.Cv) ± $(sampling.Cv_err)")
