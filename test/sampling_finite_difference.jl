let sys = DenseSystem(test_params...),
    beta = 12.34,
    dbeta = 1e-4,
    P = 256,
    num_samples = 100,
    sampling_sys = DiagonalSystem(test_params_sampling...),
    sampling_beta = 11.11

    # Only testing that it runs to completion.
    @test_nowarn SamplingFiniteDifference(sys, beta, dbeta, P, num_samples; progress_output=devnull)
    @test_nowarn SamplingFiniteDifference(sys, beta, dbeta, P, num_samples; sampling_sys=sampling_sys, sampling_beta=sampling_beta, progress_output=devnull)
end
