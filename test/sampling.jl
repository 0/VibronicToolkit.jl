using VibronicToolkit: SamplingParameters

let sys = DenseSystem(test_params...),
    beta = 12.34,
    P = 256,
    sp = SamplingParameters(sys, beta, P)

    @test sp.tau == beta/P
    @test isapprox(sp.weights/sum(sp.weights), [0.999996583, 3.41696734e-6, 1.16757213e-11])
    @test isapprox(sp.deltas, [-4.91666666e-6, -8.66666666e-6, -1.2525e-5])
    @test isapprox(sp.ds, [-1.5e-3 -1.66666666e-3 -1.75e-3; -1.33333333e-3 -1.5e-3 -1.6e-3])
end
