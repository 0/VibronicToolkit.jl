let sys = DenseSystem(test_params...)

    @test sys.freq == test_freq
    @test keys(sys.coef) == [0, 1, 2, 3, 4]
    @test sys.energy == sys.coef[0] == test_energy
    @test sys.lin == sys.coef[1] == test_lin
    @test sys.quad == sys.coef[2] == test_quad
    @test sys.coef[3] == test_cub_diag
    @test sys.coef[4] == test_quart_diag
    @test sys.coef[5] == zeros(test_M, test_M, test_M, test_M, test_M, test_S, test_S)

    @test !isdiag(sys)

    @test isdiag(diag(sys))
    @test (diag(sys)).freq == test_freq
    @test keys(diag(sys).coef) == [0, 1, 2, 3, 4]
    @test (diag(sys)).coef[0] == test_energy_diag
    @test (diag(sys)).coef[1] == test_lin_diag
    @test (diag(sys)).coef[2] == test_quad_diag
    @test (diag(sys)).coef[3] == test_cub_diag
    @test (diag(sys)).coef[4] == test_quart_diag
    @test (diag(sys)).coef[5] == zeros(test_M, test_M, test_M, test_M, test_M, test_S, test_S)

    @test simplify(sys) == DenseSystem(test_freq, [test_energy, test_lin])
    @test simplify(sys; ord=2) == DenseSystem(test_freq, [test_energy, test_lin, test_quad])
    @test simplify(sys; ord=4) == sys

    @test_throws SurfaceCouplingException DiagonalSystem(sys)
end

let sys = DenseSystem(test_params_diag...)

    @test sys.freq == test_freq
    @test keys(sys.coef) == [0, 1, 2, 3, 4]
    @test sys.energy == sys.coef[0] == test_energy_diag
    @test sys.lin == sys.coef[1] == test_lin_diag
    @test sys.quad == sys.coef[2] == test_quad_diag
    @test sys.coef[3] == test_cub_diag
    @test sys.coef[4] == test_quart_diag
    @test sys.coef[5] == zeros(test_M, test_M, test_M, test_M, test_M, test_S, test_S)

    @test isdiag(sys)

    @test diag(sys) == sys

    @test simplify(sys) == DenseSystem(test_freq, [test_energy_diag, test_lin_diag])
    @test simplify(sys; ord=2) == DenseSystem(test_freq, [test_energy_diag, test_lin_diag, test_quad_diag])
    @test simplify(sys; ord=4) == sys

    @test DiagonalSystem(sys) isa DiagonalSystem
    @test DiagonalSystem(sys) == sys
end
