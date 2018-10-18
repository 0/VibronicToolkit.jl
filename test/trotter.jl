let sys = DenseSystem(test_params...),
    beta = 12.34,
    basis_size = 10,
    P = 256,
    trotter = Trotter(sys, beta, basis_size, P)

    @test isapprox(trotter.Z, 7.78099507e-3)
    @test isapprox(trotter.E, 3.92839066e-1)
    @test isapprox(trotter.Cv, -5.05361185e-1)
end

let sys = DenseSystem(test_params...),
    trial = UniformTrialWavefunction(sys),
    beta = 12.34,
    basis_size = 10,
    P = 256,
    trotter = PigsTrotter(sys, trial, beta, basis_size, P)

    @test isapprox(trotter.Z, 2.70438291e-1)
    @test isapprox(trotter.E, 3.93991391e-1)
    @test isapprox(trotter.SvN, 9.47153746e-1)
    @test isapprox(trotter.S2, 8.27962641e-1)

    @test_throws DomainError PigsTrotter(sys, trial, beta, basis_size, P+1)
end
