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
    "--basis-size"
        metavar = "N"
        help = "single-mode basis size"
        arg_type = Int
        required = true
end
c = parse_args(ARGS, s, as_symbols=true)

sys = read(c[:conf], DenseSystem)
beta = c[:beta]
basis_size = c[:basis_size]

sos = SumOverStates(sys, beta, basis_size)
simple = Analytical(simplify(diag(sys)), beta)

println("Z/Z: $(sos.Z/simple.Z)")
println("E: $(sos.E)")
println("Cv: $(sos.Cv)")
