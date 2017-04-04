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
end
c = parsed_args = parse_args(ARGS, s, as_symbols=true)

sys = System(c[:conf])
beta = Beta(c[:beta])

analytical = Analytical(sys, beta)
simple = Analytical(simplify(sys), beta)

println("Z/Z: $(analytical.Z/simple.Z)")
println("E: $(analytical.E)")
println("Cv: $(analytical.Cv)")
