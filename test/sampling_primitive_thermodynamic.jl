for (num_boot, num_seqs) in [(nothing, nothing), (10, nothing), (10, 4)]
    let sys = DenseSystem(test_params...),
        beta = 12.34,
        P = 256,
        num_samples = 100,
        sm = make_sampling_method(num_samples, num_boot, num_seqs),
        sampling_sys = DiagonalSystem(test_params_sampling...)

        # Only testing that it runs to completion.
        @test_nowarn SamplingPrimitiveThermodynamic(sys, beta, P, sm; progress_output=devnull)
        @test_nowarn SamplingPrimitiveThermodynamic(sys, beta, P, sm; sampling_sys=sampling_sys, progress_output=devnull)
    end
end
