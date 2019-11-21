let sys = DenseSystem(test_params_flat...),
    max_iter = 100,
    decomp = IterativeDecomposition(sys, max_iter),
    sys_decomp = DiagonalSystem(decomp)

    @test decomp.degens == [1]

    # Remove overall phase.
    vs_decomp = decomp.vs[1,1] > 0 ? decomp.vs : -decomp.vs
    vs_expected = reshape([0.5833308951, 0.5773294744, 0.5713280536], 3, 1)
    @test isapprox(vs_decomp, vs_expected)

    expected_energy = reshape([-2.880208317], 1, 1)
    expected_freq = reshape([1.0, 2.0], test_M, 1)
    expected_lin = reshape([0.0149573419, 0.0179571258], test_M, 1, 1)
    expected_quad = reshape([8.97856291e-3 1.04784549e-2; 1.04784549e-2 1.19783468e-2], test_M, test_M, 1, 1)
    sys_expected = DiagonalSystem(expected_freq, [expected_energy, expected_lin, expected_quad])
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
