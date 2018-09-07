let sys = DenseSystem(test_params...),
    beta = 12.34,
    dbeta = 1e-4,
    P = 256,
    num_samples = 100,
    sampling = SamplingFiniteDifference(sys, beta, dbeta, P, num_samples; progress_output=devnull)

    # Only testing that it runs to completion.
end
