#!/usr/bin/env julia

using VibronicToolkit

using ArgParse
using PyPlot

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
    "--conf"
        metavar = "FILE"
        help = "path to config file"
        required = true
    "--sampling-conf"
        metavar = "FILE"
        help = "path to sampling config file"
    "--sampling-beta"
        metavar = "T"
        help = "sampling reciprocal temperature"
        arg_type = Float64
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
    "--contour"
        help = "draw contours instead of smooth gradients"
        action = :store_true
end
c = parse_args(ARGS, s, as_symbols=true)

sys = read(c[:conf], DenseSystem)
if c[:sampling_conf] !== nothing
    sampling_sys = read(c[:sampling_conf], DiagonalSystem)
    sampling_beta = c[:sampling_beta]
    P = c[:num_links]
end
extent = Tuple(parse.(Float64, split(c[:extent], ",")))
out_path = c[:out_path]

# Draw PES.
gp = GroundPes(sys, extent)
pes = gp.pes
if c[:contour]
    cf = contourf(pes, origin="lower", extent=extent, cmap="magma_r")
    contour(pes, origin="lower", extent=extent, colors="0.5", linestyles="solid", linewidths=0.2)
else
    cf = imshow(pes, origin="lower", aspect="auto", extent=extent, cmap="magma_r")
end
colorbar(cf)

if c[:sampling_conf] !== nothing
    # Draw sampling distribution.
    ds = sampling_sys.ds
    ws = weights(sampling_sys, sampling_beta)
    colors = map(plt.matplotlib.colors.hsv_to_rgb, [(0.4, w, 0.8) for w in ws])
    scatter(ds[1,:], ds[2,:], marker="+", c=colors)
    for s in 1:size(ds, 2)
        sigma_1 = path_mean_std(sampling_sys, sampling_beta, P, s, 1)
        sigma_2 = path_mean_std(sampling_sys, sampling_beta, P, s, 2)
        ellipse_f = plt.matplotlib.patches.Ellipse
        ellipse = ellipse_f(xy=(ds[1,s], ds[2,s]), width=2sigma_1, height=2sigma_2, edgecolor=colors[s], fc="None", linestyle="dotted")
        gca().add_patch(ellipse)
        offset = (extent[2] - extent[1]) / 32
        text(ds[1,s]+offset, ds[2,s], "$(s)", ha="center", va="center", color="0.5")
    end
end

xlabel(L"q_1")
ylabel(L"q_2")

savefig(out_path, bbox_inches="tight")
