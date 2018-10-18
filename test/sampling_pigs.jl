let sys = DenseSystem(test_params...),
    beta = 12.34,
    P = 256,
    num_samples = 100,
    sampling_sys = DiagonalSystem(test_params_sampling...)

    # Only testing that it runs to completion.
    @test_nowarn PigsSampling(sys, beta, P, num_samples; progress_output=devnull)
    @test_nowarn PigsSampling(sys, beta, P, num_samples; sampling_sys=sampling_sys, progress_output=devnull)
end
