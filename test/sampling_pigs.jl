for num_boot in [nothing, 10]
    let sys = DenseSystem(test_params...),
        trial = UniformTrialWavefunction(sys),
        beta = 12.34,
        P = 256,
        num_samples = 100,
        sm = make_sampling_method(num_samples, num_boot),
        sampling_sys = DiagonalSystem(test_params_sampling...),
        sampling_trial = UniformTrialWavefunction(sampling_sys)

        # Only testing that it runs to completion.
        @test_nowarn PigsSampling(sys, trial, beta, P, sm; progress_output=devnull)
        @test_nowarn PigsSampling(sys, trial, beta, P, sm; sampling_sys=sampling_sys, sampling_trial=sampling_trial, progress_output=devnull)
    end
end
