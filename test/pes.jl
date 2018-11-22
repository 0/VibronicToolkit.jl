let sys = DenseSystem(test_params...),
    beta = 12.34,
    P = 256,
    sampling_sys = DiagonalSystem(test_params_sampling...),
    gp = GroundPes(sys, (-1.0, 2.0, -3.0, 4.0); lengths=(2, 3))

    @test isapprox(gp.pes, [1.34545191e1 1.62876293e1; -9.81542994e-1 2.81086523; 2.40243700e1 2.83118307e1])

    @test isapprox(path_mean_std(sampling_sys, beta, P, 1, 1), 0.201370345)
    @test isapprox(path_mean_std(sampling_sys, beta, P, test_S_sampling, test_M), 0.108103579)
end
