let sys = DiagonalSystem(test_params_diag...)

    @test sys.energy == test_energy_diag
    @test sys.freq == test_freq
    @test sys.lin == test_lin_diag
    @test sys.quad == test_quad_diag

    @test isdiag(sys)

    @test diag(sys) == sys

    @test simplify(sys) == DiagonalSystem(test_energy_diag, test_freq, test_lin_diag, zero(test_quad_diag))

    @test DenseSystem(sys) isa DenseSystem
    @test DenseSystem(sys) == sys
end

@test_throws DomainError DiagonalSystem(test_freq, test_energy_diag, test_lin_diag, test_quad_diag)
@test_throws SurfaceCouplingException DiagonalSystem(test_params...)
