let sys = DiagonalSystem(test_params_diag...)

    @test sys.freq == test_freq
    @test keys(sys.coef) == [0, 1, 2]
    @test sys.energy == sys.coef[0] == test_energy_diag
    @test sys.lin == sys.coef[1] == test_lin_diag
    @test sys.quad == sys.coef[2] == test_quad_diag
    @test sys.coef[3] == zeros(test_M, test_M, test_M, test_S, test_S)

    @test isapprox(sys.deltas, [-4.91666666e-6, -8.66666666e-6, -1.2525e-5])
    @test isapprox(sys.ds, [-1.5e-3 -1.66666666e-3 -1.75e-3; -1.33333333e-3 -1.5e-3 -1.6e-3])

    @test isdiag(sys)

    @test diag(sys) == sys

    @test simplify(sys) == DiagonalSystem(test_freq, [test_energy_diag, test_lin_diag])
    @test simplify(sys; ord=2) == sys

    @test DenseSystem(sys) isa DenseSystem
    @test DenseSystem(sys) == sys
end

@test_throws SurfaceCouplingException DiagonalSystem(test_params...)
