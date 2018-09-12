#!/usr/bin/env julia

using VibronicToolkit

using DelimitedFiles

using ArgParse

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
    "--conf"
        metavar = "FILE"
        help = "path to config file"
        required = true
    "--max-iter"
        metavar = "ITER"
        help = "maximum number of iterations"
        arg_type = Int
        default = 100
    "--out-conf"
        metavar = "FILE"
        help = "path to system configuration output file"
        required = true
    "--out-vs"
        metavar = "FILE"
        help = "path to transformation vector output file"
end
c = parse_args(ARGS, s, as_symbols=true)

sys = read(c[:conf], DenseSystem)
max_iter = c[:max_iter]
out_conf = c[:out_conf]
out_vs = c[:out_vs]

decomp = IterativeDecomposition(sys, max_iter)

write(out_conf, DiagonalSystem(decomp))
if out_vs !== nothing
    writedlm(out_vs, permutedims(decomp.vs))
end
