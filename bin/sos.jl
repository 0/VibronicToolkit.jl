#!/usr/bin/env julia

using VibronicToolkit

using ArgParse

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
    "--pigs"
        help = "use PIGS instead of finite temperature"
        action = :store_true
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

if c[:pigs]
    trial = UniformTrialWavefunction(sys)
    sos = PigsSumOverStates(sys, trial, beta, basis_size)

    println("Z: $(sos.Z)")
    println("E: $(sos.E)")
    println("E0_exact: $(sos.E0_exact)")
else
    sos = SumOverStates(sys, beta, basis_size)

    println("Z: $(sos.Z)")
    println("E: $(sos.E)")
    println("Cv: $(sos.Cv)")
end
