let sys = DiagonalSystem(test_params_diag...),
    beta = 12.34,
    analytical = Analytical(sys, beta)

    @test isapprox(analytical.Z, 6.92896484e-9)
    @test isapprox(analytical.E, 1.52249605e0)
    @test isapprox(analytical.Cv, 5.35769920e-4)
end

let sys = DiagonalSystem(test_params_diag...),
    beta = 12.34,
    analytical = PigsAnalytical(sys, beta)

    @test isapprox(analytical.Z, 8.70719401e-8)
    @test isapprox(analytical.E, 1.52249605e0)
end
