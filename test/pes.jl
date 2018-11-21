let sys = DenseSystem(test_params...),
    beta = 12.34,
    P = 256,
    sampling_sys = DiagonalSystem(test_params_sampling...),
    gp = GroundPes(sys, (-1.0, 2.0, -3.0, 4.0); lengths=(2, 3))

    @test isapprox(gp.pes, [13.2683190 16.2872293; -0.98152163 2.77276916; 23.9078621 27.0343054])

    @test isapprox(path_mean_std(sampling_sys, beta, P, 1, 1), 0.201370345)
    @test isapprox(path_mean_std(sampling_sys, beta, P, test_S_sampling, test_M), 0.108103579)
end
