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
    "--sampling-conf"
        metavar = "FILE"
        help = "path to sampling config file"
    "--sampling-beta"
        metavar = "T"
        help = "sampling reciprocal temperature"
        arg_type = Float64
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
    # Changing the sampling system changes the normalization for Z, so we
    # renormalize it to obtain the same result as when using the default
    # sampling system.
    simple_reg = Analytical(simplify(diag(sys)), beta)
    simple_samp = Analytical(sampling_sys, beta)
    Z_renorm = simple_samp.Z / simple_reg.Z
else
    sampling_sys = nothing
    Z_renorm = 1.0
end
sampling_beta = c[:sampling_beta]

if dbeta !== nothing
    sampling = SamplingFiniteDifference(sys, beta, dbeta, P, num_samples; sampling_sys=sampling_sys, sampling_beta=sampling_beta)
else
    sampling = SamplingPrimitiveThermodynamic(sys, beta, P, num_samples; sampling_sys=sampling_sys, sampling_beta=sampling_beta)
end

println("Z/Z: $(sampling.Z * Z_renorm) ± $(sampling.Z_err * Z_renorm)")
println("E: $(sampling.E) ± $(sampling.E_err)")
println("Cv: $(sampling.Cv) ± $(sampling.Cv_err)")
