let sys = DenseSystem(test_params...)

    @test sys.energy == test_energy
    @test sys.freq == test_freq
    @test sys.lin == test_lin
    @test sys.quad == test_quad

    @test !isdiag(sys)

    @test isdiag(diag(sys))
    @test (diag(sys)).energy == test_energy_diag
    @test (diag(sys)).freq == test_freq
    @test (diag(sys)).lin == test_lin_diag
    @test (diag(sys)).quad == test_quad_diag

    @test simplify(sys) == DenseSystem(test_energy, test_freq, test_lin, zero(test_quad))

    @test_throws SurfaceCouplingException DiagonalSystem(sys)
end

let sys = DenseSystem(test_params_diag...)

    @test sys.energy == test_energy_diag
    @test sys.freq == test_freq
    @test sys.lin == test_lin_diag
    @test sys.quad == test_quad_diag

    @test isdiag(sys)

    @test diag(sys) == sys

    @test simplify(sys) == DenseSystem(test_energy_diag, test_freq, test_lin_diag, zero(test_quad_diag))

    @test DiagonalSystem(sys) isa DiagonalSystem
    @test DiagonalSystem(sys) == sys
end

@test_throws DomainError DenseSystem(test_freq, test_energy, test_lin, test_quad)
