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
    "--uncoupled"
        help = "remove inter-surface and quadratic coupling from system"
        action = :store_true
end
c = parse_args(ARGS, s, as_symbols=true)

sys = System(c[:conf])
beta = c[:beta]

simple = Analytical(simplify(sys), beta)

if c[:uncoupled]
    println("Z: $(simple.Z)")
    println("E: $(simple.E)")
    println("Cv: $(simple.Cv)")
else
    analytical = Analytical(sys, beta)

    println("Z/Z: $(analytical.Z/simple.Z)")
    println("E: $(analytical.E)")
    println("Cv: $(analytical.Cv)")
end
