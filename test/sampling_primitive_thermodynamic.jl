let sys = DenseSystem(test_params...),
    beta = 12.34,
    P = 256,
    num_samples = 100,
    sampling = SamplingPrimitiveThermodynamic(sys, beta, P, num_samples; progress_output=devnull)

    # Only testing that it runs to completion.
end
