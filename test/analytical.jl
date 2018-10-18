let sys = DiagonalSystem(test_params_diag...),
    beta = 12.34,
    analytical = Analytical(sys, beta)

    @test isapprox(analytical.Z, 6.92896484e-9)
    @test isapprox(analytical.E, 1.52249605e0)
    @test isapprox(analytical.Cv, 5.35769920e-4)
end

let sys = DiagonalSystem(test_params_diag...),
    trial = UniformTrialWavefunction(sys),
    beta = 12.34,
    analytical = PigsAnalytical(sys, trial, beta)

    @test isapprox(analytical.Z, 8.70719401e-8)
    @test isapprox(analytical.E, 1.52249605e0)
    @test isapprox(analytical.SvN, 4.7216113e-11)
    @test isapprox(analytical.S2, 3.3582026e-12)
end

let sys = DiagonalSystem(test_params_diag...),
    trial = UniformTrialWavefunction(sys),
    beta = 1e-5
    analytical = PigsAnalytical(sys, trial, beta)

    @test isapprox(analytical.Z, 577814.838)
    @test isapprox(analytical.E, 99999.0360)
    @test isapprox(analytical.SvN, 7.87900924e-2)
    @test isapprox(analytical.S2, 3.01336884e-2)
end
