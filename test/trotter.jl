let sys = DenseSystem(test_params...),
    beta = 12.34,
    basis_size = 10,
    P = 256,
    trotter = Trotter(sys, beta, basis_size, P)

    @test isapprox(trotter.Z, 7.52871676e-3)
    @test isapprox(trotter.E, 3.95509128e-1)
    @test isapprox(trotter.Cv, -5.07311833e-1)
end

let sys = DenseSystem(test_params...),
    trial = UniformTrialWavefunction(sys),
    beta = 12.34,
    basis_size = 10,
    P = 256,
    trotter = PigsTrotter(sys, trial, beta, basis_size, P; splitting=p2_U)

    @test isapprox(trotter.Z, 2.71333560e-1)
    @test isapprox(trotter.E, 3.83343408e-1)
    @test isapprox(trotter.SvN, 9.47734102e-1)
    @test isapprox(trotter.S2, 8.28872984e-1)

    @test_throws DomainError PigsTrotter(sys, trial, beta, basis_size, P+1)
end
