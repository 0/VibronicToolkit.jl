using VibronicToolkit: Basis, trial_mode, trial_spatial

let sys = DenseSystem(test_params_trivial...),
    size = 3,
    basis = Basis(sys, size),
    surface_coefs = [5.0, 2.0, 0.5],
    trial = UniformTrialWavefunction(sys, surface_coefs)

    @test isapprox(trial_spatial(trial, [9.0, 11.0]), surface_coefs)

    @test isapprox(trial_mode(trial, basis), [17.7245385, 0.0, 12.5331413, 0.0,
                                              0.0, 0.0, 12.5331413, 0.0,
                                              8.86226925,
                                              7.08981540, 0.0, 5.01325654, 0.0,
                                              0.0, 0.0, 5.01325654, 0.0,
                                              3.54490770,
                                              1.77245385, 0.0, 1.25331414, 0.0,
                                              0.0, 0.0, 1.25331414, 0.0,
                                              0.886226925])
end
