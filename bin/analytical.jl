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

sys = read(c[:conf], DenseSystem)
beta = c[:beta]

if c[:uncoupled]
    sys = simplify(diag(sys))
else
    try
        global sys = DiagonalSystem(sys)
    catch ex
        ex isa SurfaceCouplingException || rethrow()
        error("Analytical solution only applies to uncoupled systems")
    end
end

analytical = Analytical(sys, beta)

println("Z: $(analytical.Z)")
println("E: $(analytical.E)")
println("Cv: $(analytical.Cv)")
