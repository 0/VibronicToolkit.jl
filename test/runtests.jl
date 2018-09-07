using VibronicToolkit
using Test

test_S = 3
test_M = 2
test_energy = zeros(test_S, test_S)
test_energy_diag = zeros(test_S, test_S)
test_freq = zeros(test_M, test_S)
# The JSON format doesn't support different frequencies between surfaces.
test_freq_flat = zeros(test_M, test_S)
test_lin = zeros(test_M, test_S, test_S)
test_lin_diag = zeros(test_M, test_S, test_S)
test_quad = zeros(test_M, test_M, test_S, test_S)
test_quad_diag = zeros(test_M, test_M, test_S, test_S)
for s1 in 1:test_S
    for m1 in 1:test_M
        test_freq[m1, s1] = s1+m1
        test_freq_flat[m1, s1] = m1
    end
    test_energy_diag[s1, s1] = 0.01*(s1+s1) - 1.0
    for m1 in 1:test_M
        test_lin_diag[m1, s1, s1] = 0.001*(s1+s1+m1)
        for m2 in 1:test_M
            test_quad_diag[m2, m1, s1, s1] = 0.001*(s1+s1+m1+m2)
        end
    end
    for s2 in 1:test_S
        test_energy[s2, s1] = 0.01*(s1+s2) - 1.0
        for m1 in 1:test_M
            test_lin[m1, s2, s1] = 0.001*(s1+s2+m1)
            for m2 in 1:test_M
                test_quad[m2, m1, s2, s1] = 0.001*(s1+s2+m1+m2)
            end
        end
    end
end
test_params = (test_energy, test_freq, test_lin, test_quad)
test_params_diag = (test_energy_diag, test_freq, test_lin_diag, test_quad_diag)
test_params_flat = (test_energy, test_freq_flat, test_lin, test_quad)
test_params_flat_diag = (test_energy_diag, test_freq_flat, test_lin_diag, test_quad_diag)
test_params_trivial = (ones(test_S, test_S),
                       ones(test_M, test_S),
                       zeros(test_M, test_S, test_S),
                       zeros(test_M, test_M, test_S, test_S))

@testset "VibronicToolkit" begin
    @testset "utilities.jl" begin include("utilities.jl") end

    @testset "system.jl" begin include("system.jl") end
    @testset "operators.jl" begin include("operators.jl") end

    @testset "jackknife.jl" begin include("jackknife.jl") end

    @testset "analytical.jl" begin include("analytical.jl") end
    @testset "sos.jl" begin include("sos.jl") end
    @testset "trotter.jl" begin include("trotter.jl") end

    @testset "sampling.jl" begin include("sampling.jl") end
    @testset "sampling_finite_difference.jl" begin include("sampling_finite_difference.jl") end
    @testset "sampling_primitive_thermodynamic.jl" begin include("sampling_primitive_thermodynamic.jl") end
end
