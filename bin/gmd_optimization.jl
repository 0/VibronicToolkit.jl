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
    "--num-boot"
        metavar = "B"
        help = "number of bootstrap samples"
        arg_type = Int
    "--num-seqs"
        metavar = "S"
        help = "number of randomized quasi-random sequences"
        arg_type = Int
    "--num-iter"
        metavar = "N"
        help = "number of SPSA iterations"
        arg_type = Int
        required = true
    "--spsa-a"
        metavar = "A"
        help = "SPSA gain parameter"
        arg_type = Float64
        required = true
    "--loss-bound-lower"
        metavar = "LL"
        help = "lower bound on loss function"
        arg_type = Float64
        required = true
    "--loss-bound-upper"
        metavar = "LU"
        help = "upper bound on loss function"
        arg_type = Float64
        required = true
    "--repulsion-stiffness"
        metavar = "S"
        help = "stiffness in repulsion term"
        arg_type = Float64
        required = true
    "--repulsion-coefficient"
        metavar = "R"
        help = "coefficient of repulsion term"
        arg_type = Float64
        required = true
    "--weight-imbalance-coefficient"
        metavar = "W"
        help = "coefficient of weight imbalance term"
        arg_type = Float64
        required = true
    "--num-surfaces"
        metavar = "S"
        help = "number of GMD surfaces"
        arg_type = Int
    "--start-conf"
        metavar = "FILE"
        help = "path to starting config file"
    "--out-data"
        metavar = "FILE"
        help = "path to data output file"
    "--out-conf"
        metavar = "FILE"
        help = "path to optimized config file"
        required = true
    "--quiet"
        help = "suppress progress meter"
        action = :store_true
end
c = parse_args(ARGS, s, as_symbols=true)

sys = read(c[:conf], DenseSystem)
beta = c[:beta]
P = c[:num_links]
num_samples = c[:num_samples]
num_boot = c[:num_boot]
num_seqs = c[:num_seqs]
num_iter = c[:num_iter]
spsa_a = c[:spsa_a]
num_surfaces = c[:num_surfaces]
if !isnothing(c[:out_data])
    out_data_file = open(c[:out_data], "w")
else
    out_data_file = devnull
end
out_conf = c[:out_conf]
if c[:quiet]
    progress_output = devnull
else
    progress_output = stderr
end

sm = make_sampling_method(num_samples, num_boot, num_seqs)
gop = GmdOptimizationParameters{P}(beta, sm, c[:loss_bound_lower],
                                   c[:loss_bound_upper],
                                   c[:repulsion_stiffness],
                                   c[:repulsion_coefficient],
                                   c[:weight_imbalance_coefficient])

if !isnothing(c[:start_conf])
    start_sys = read(c[:start_conf], DiagonalSystem)
    opt = GmdOptimization(sys, gop, num_iter, spsa_a, start_sys;
                          progress_output=progress_output,
                          data_output=out_data_file)
elseif !isnothing(num_surfaces)
    opt = GmdOptimizationDeformation(sys, gop, num_iter, spsa_a, num_surfaces;
                                     progress_output=progress_output,
                                     data_output=out_data_file)
else
    error("Start system or number of surfaces must be given")
end

write(out_conf, DiagonalSystem(opt))
