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
    "--basis-size"
        metavar = "N"
        help = "single-mode basis size"
        arg_type = Int
        required = true
    "--num-links"
        metavar = "P"
        help = "number of Trotter links"
        arg_type = Int
        required = true
end
c = parse_args(ARGS, s, as_symbols=true)

sys = System(c[:conf])
beta = c[:beta]
basis_size = c[:basis_size]
P = c[:num_links]

trotter = Trotter(sys, beta, basis_size, P)
simple = Analytical(simplify(sys), beta)

println("Z/Z: $(trotter.Z/simple.Z)")
println("E: $(trotter.E)")
println("Cv: $(trotter.Cv)")
