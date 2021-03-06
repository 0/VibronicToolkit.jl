let sys = DenseSystem(test_params...),
    beta = 12.34,
    basis_size = 10,
    sos = SumOverStates(sys, beta, basis_size)

    @test isapprox(sos.Z, 7.74822371e-3)
    @test isapprox(sos.E, 3.93864802e-1)
    @test isapprox(sos.Cv, 1.04927608e-10)
end
