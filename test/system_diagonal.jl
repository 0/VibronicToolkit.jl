let sys = DiagonalSystem(test_params_diag...)

    @test sys.energy == test_energy_diag
    @test sys.freq == test_freq
    @test sys.lin == test_lin_diag
    @test sys.quad == test_quad_diag

    @test isapprox(sys.deltas, [-4.91666666e-6, -8.66666666e-6, -1.2525e-5])
    @test isapprox(sys.ds, [-1.5e-3 -1.66666666e-3 -1.75e-3; -1.33333333e-3 -1.5e-3 -1.6e-3])

    @test isdiag(sys)

    @test diag(sys) == sys

    @test simplify(sys) == DiagonalSystem(test_energy_diag, test_freq, test_lin_diag, zero(test_quad_diag))

    @test DenseSystem(sys) isa DenseSystem
    @test DenseSystem(sys) == sys
end

@test_throws DomainError DiagonalSystem(test_freq, test_energy_diag, test_lin_diag, test_quad_diag)
@test_throws SurfaceCouplingException DiagonalSystem(test_params...)
