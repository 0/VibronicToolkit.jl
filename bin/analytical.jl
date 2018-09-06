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

sys = DenseSystem(c[:conf])
beta = c[:beta]

simple = Analytical(simplify(diag(sys)), beta)

if c[:uncoupled]
    println("Z: $(simple.Z)")
    println("E: $(simple.E)")
    println("Cv: $(simple.Cv)")
else
    try
        global analytical = Analytical(DiagonalSystem(sys), beta)
    catch ex
        ex isa SurfaceCouplingException || rethrow()
        error("Analytical solution only applies to uncoupled systems")
    end

    println("Z/Z: $(analytical.Z/simple.Z)")
    println("E: $(analytical.E)")
    println("Cv: $(analytical.Cv)")
end
