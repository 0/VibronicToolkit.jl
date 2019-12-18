let sys = DenseSystem(test_params...),
    beta = 12.34,
    P = 4,
    num_samples = 100,
    sm = make_sampling_method(num_samples, nothing, nothing),
    gop = GmdOptimizationParameters{P}(beta, sm, 1.0, 2.0, 3.0, 1.0, 1.5),
    num_iter = 2,
    spsa_a = 1e-2,
    start_sys = DiagonalSystem(test_params_sampling...),
    num_surfaces = test_S + 1

    # Only testing that it runs to completion.
    @test_nowarn GmdOptimization(sys, gop, num_iter, spsa_a, start_sys; progress_output=devnull)
    @test_nowarn GmdOptimizationDeformation(sys, gop, num_iter, spsa_a, num_surfaces; progress_output=devnull)
end
