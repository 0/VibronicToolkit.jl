#!/usr/bin/env julia

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "../src"))
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
c = parsed_args = parse_args(ARGS, s, as_symbols=true)

sys = System(c[:conf])
beta = c[:beta]
P = c[:num_links]
num_samples = c[:num_samples]

sampling = Sampling(sys, beta, P, num_samples)

println("Z/Z: $(sampling.Z) ± $(sampling.Z_err)")
