#!/usr/bin/env julia

using VibronicToolkit

using ArgParse
using PyPlot

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
    "--sampling-conf"
        metavar = "FILE"
        help = "path to sampling config file"
    "--num-links"
        metavar = "P"
        help = "number of Trotter links"
        arg_type = Int
    "--extent"
        metavar = "q1_min,q1_max,q2_min,q2_max"
        help = "comma-separated extent values"
        required = true
    "--out-path"
        metavar = "FILE"
        help = "path to output file"
        required = true
    "--out-path-exact"
        metavar = "FILE"
        help = "path to output file for exact density (PIGS only)"
    "--sum-diagonal"
        help = "only generate overall density"
        action = :store_true
    "--contour"
        help = "draw contours instead of smooth gradients"
        action = :store_true
    "--quiet"
        help = "suppress progress meter"
        action = :store_true
end
c = parse_args(ARGS, s, as_symbols=true)

sys = read(c[:conf], DenseSystem)
beta = c[:beta]
basis_size = c[:basis_size]
if c[:sampling_conf] !== nothing
    sampling_sys = read(c[:sampling_conf], DiagonalSystem)
    P = c[:num_links]
end
extent = Tuple(parse.(Float64, split(c[:extent], ",")))
out_path = c[:out_path]
if c[:quiet]
    progress_output = devnull
else
    progress_output = stderr
end

if c[:pigs]
    trial = UniformTrialWavefunction(sys)
    density_all, density_all_exact = diagonal_density_pigs(sys, trial, beta, basis_size, extent; progress_output=progress_output)
else
    density_all = diagonal_density(sys, beta, basis_size, extent; progress_output=progress_output)
end

densities = [density_all]
out_paths = [out_path]
if c[:pigs] && c[:out_path_exact] !== nothing
    push!(densities, density_all_exact)
    push!(out_paths, c[:out_path_exact])
end

for (density_all, out_path) in zip(densities, out_paths)
    if c[:sum_diagonal]
        density_sum = zero(density_all[:, :, 1:1, 1:1])
        for s in 1:size(density_all, 4)
            density_sum[:, :, 1, 1] += density_all[:, :, s, s]
        end
        density_all = density_sum
    end
    vmin, vmax = minimum(density_all), maximum(density_all)

    if c[:contour]
        levels = range(vmin; stop=vmax, length=13)
    end

    clf()
    fig, axes = subplots(size(density_all, 3), size(density_all, 4); sharex=true, sharey=true, squeeze=false)

    for s1 in 1:size(density_all, 4)
        for s2 in 1:size(density_all, 3)
            density = density_all[:, :, s2, s1]
            ax = axes[s2, s1]

            # Draw density.
            if c[:contour]
                global cf = ax["contourf"](density, origin="lower", extent=extent, levels=levels, cmap="magma")
                if findfirst(x -> minimum(density) <= x, levels) < findlast(x -> x <= maximum(density), levels)
                    # Only show the contours if any appear in the current
                    # density landscape.
                    ax["contour"](density, origin="lower", extent=extent, levels=levels, colors="0.5", linestyles="solid", linewidths=0.2)
                end
            else
                global cf = ax["imshow"](density, origin="lower", aspect="auto", extent=extent, vmin=vmin, vmax=vmax, cmap="magma")
            end

            if c[:sampling_conf] !== nothing
                # Draw sampling distribution.
                ds = sampling_sys.ds
                ws = weights(sampling_sys, beta)
                colors = map(plt["matplotlib"]["colors"]["hsv_to_rgb"], [(0.4, w, 0.8) for w in ws])
                ax["scatter"](ds[1,:], ds[2,:], marker="+", c=colors)
                for s in 1:size(ds, 2)
                    sigma_1 = path_mean_std(sampling_sys, beta, P, s, 1)
                    sigma_2 = path_mean_std(sampling_sys, beta, P, s, 2)
                    ellipse_f = plt["matplotlib"]["patches"]["Ellipse"]
                    ellipse = ellipse_f(xy=(ds[1,s], ds[2,s]), width=2sigma_1, height=2sigma_2, edgecolor=colors[s], fc="None", linestyle="dotted")
                    ax["add_patch"](ellipse)
                    offset = (extent[2] - extent[1]) * size(density_all, 4) / 32
                    ax["text"](ds[1,s]+offset, ds[2,s], "$(s)", ha="center", va="center", color="0.5")
                end
            end
        end
    end

    for s1 in 1:size(density_all, 4)
        axes[end, s1]["set_xlabel"](L"q_1")
    end

    for s2 in 1:size(density_all, 3)
        axes[s2, 1]["set_ylabel"](L"q_2")
    end

    fig["colorbar"](cf; ax=vec(axes), use_gridspec=true)

    savefig(out_path, bbox_inches="tight")
end
