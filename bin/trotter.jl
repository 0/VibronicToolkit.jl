#!/usr/bin/env julia

using VibronicToolkit

using ArgParse

function ArgParse.parse_item(::Type{Splitting}, x::AbstractString)
    if x == "h0_V"
        h0_V
    elseif x == "p2_U"
        p2_U
    else
        error("Invalid splitting: $(x)")
    end
end

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
    "--num-links"
        metavar = "P"
        help = "number of Trotter links"
        arg_type = Int
        required = true
    "--splitting"
        metavar = "S"
        help = "choice of operator splitting (h0_V, p2_U)"
        arg_type = Splitting
        default = h0_V
end
c = parse_args(ARGS, s, as_symbols=true)

sys = read(c[:conf], DenseSystem)
beta = c[:beta]
basis_size = c[:basis_size]
P = c[:num_links]
splitting = c[:splitting]

if c[:pigs]
    trial = UniformTrialWavefunction(sys)
    trotter = PigsTrotter(sys, trial, beta, basis_size, P; splitting=splitting)

    println("Z: $(trotter.Z)")
    println("E: $(trotter.E)")
    println("SvN: $(trotter.SvN)")
    println("S2: $(trotter.S2)")
else
    trotter = Trotter(sys, beta, basis_size, P; splitting=splitting)

    println("Z: $(trotter.Z)")
    println("E: $(trotter.E)")
    println("Cv: $(trotter.Cv)")
end
