let sys = DenseSystem(test_params_flat...),
    max_iter = 100,
    decomp = IterativeDecomposition(sys, max_iter),
    sys_decomp = DiagonalSystem(decomp)

    @test decomp.degens == [1]

    # Remove overall phase.
    vs_decomp = decomp.vs[1,1] > 0 ? decomp.vs : -decomp.vs
    vs_expected = reshape([0.5833302853, 0.5773294386, 0.5713287124], 3, 1)
    @test isapprox(vs_decomp, vs_expected)

    expected_energy = reshape([-2.880208317], 1, 1)
    expected_freq = reshape([1.0, 2.0], test_M, 1)
    expected_lin = reshape([0.0149573465, 0.0179571305], test_M, 1, 1)
    expected_quad = reshape([8.97856524e-3, 1.04784572e-2, 1.04784572e-2, 1.19783492e-2],
                            test_M, test_M, 1, 1)
    expected_cub = reshape([6.97228455e-4, 7.97228455e-4, 7.97228455e-4, 8.97228455e-4,
                            7.97228455e-4, 8.97228455e-4, 8.97228455e-4, 9.97228455e-4],
                           test_M, test_M, test_M, 1, 1)
    expected_quart = reshape([7.97228455e-4, 8.97228455e-4, 8.97228455e-4, 9.97228455e-4,
                              8.97228455e-4, 9.97228455e-4, 9.97228455e-4, 1.09722846e-3,
                              8.97228455e-4, 9.97228455e-4, 9.97228455e-4, 1.09722846e-3,
                              9.97228455e-4, 1.09722846e-3, 1.09722846e-3, 1.19722846e-3],
                             test_M, test_M, test_M, test_M, 1, 1)
    sys_expected = DiagonalSystem(expected_freq,
                                  [expected_energy, expected_lin,
                                   expected_quad, expected_cub,
                                   expected_quart])
    @test isapprox(sys_decomp, sys_expected)
end

let sys = DenseSystem(test_params_sym...),
    max_iter = 100,
    decomp = IterativeDecomposition(sys, max_iter)

    @test decomp.degens == [1, 2]
end

let sys = DenseSystem(test_params_sym3...),
    max_iter = 100

    @test_throws DegeneracyException IterativeDecomposition(sys, max_iter)
end
