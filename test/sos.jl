let sys = DenseSystem(test_params...),
    beta = 12.34,
    basis_size = 10,
    sos = SumOverStates(sys, beta, basis_size)

    @test isapprox(sos.Z, 7.49696629e-3)
    @test isapprox(sos.E, 3.96536210e-1)
    @test isapprox(sos.Cv, 1.00260246e-10)
end

let sys = DenseSystem(test_params...),
    trial = UniformTrialWavefunction(sys),
    beta = 12.34,
    basis_size = 10,
    sos = PigsSumOverStates(sys, trial, beta, basis_size)

    @test isapprox(sos.Z, 2.59946834e-1)
    @test isapprox(sos.E, 3.96536210e-1)
    @test isapprox(sos.SvN, 9.46722467e-1)
    @test isapprox(sos.S2, 8.27305940e-1)
    @test isapprox(sos.E0_exact, 3.96536210e-1)
    @test isapprox(sos.SvN_exact, 9.46722453e-1)
    @test isapprox(sos.S2_exact, 8.27305915e-1)
end
