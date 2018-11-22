let sys = DenseSystem(test_params...),
    beta = 12.34,
    basis_size = 10,
    dd = DiagonalDensity(sys, beta, basis_size, (-1.0, 2.0, -3.0, 4.0); lengths=(2, 3), progress_output=devnull)

    @test isapprox(dd.density[:, :, 1, 1], [8.11161311e-6 4.25498516e-7; 5.40791485e-2 2.53327309e-3; 7.24246880e-9 2.54194630e-10])
    @test isapprox(dd.density[:, :, 2, 3], [2.73374202e-6 1.43505846e-7; 1.82083624e-2 8.47432180e-4; 2.40797685e-9 8.48749924e-11])
end

let sys = DenseSystem(test_params...),
    trial = UniformTrialWavefunction(sys),
    beta = 12.34,
    basis_size = 10,
    dd = PigsDiagonalDensity(sys, trial, beta, basis_size, (-1.0, 2.0, -3.0, 4.0); lengths=(2, 3), progress_output=devnull)

    @test isapprox(dd.density[:, :, 1, 1], [8.11161291e-6 4.25498503e-7; 5.40791472e-2 2.53327301e-3; 7.24246861e-9 2.54194621e-10])
    @test isapprox(dd.density[:, :, 2, 3], [2.73374211e-6 1.43505850e-7; 1.82083630e-2 8.47432207e-4; 2.40797694e-9 8.48749951e-11])
    @test isapprox(dd.density_exact[:, :, 1, 1], [8.11161311e-6 4.25498516e-7; 5.40791485e-2 2.53327309e-3; 7.24246880e-9 2.54194630e-10])
    @test isapprox(dd.density_exact[:, :, 2, 3], [2.73374202e-6 1.43505846e-7; 1.82083624e-2 8.47432180e-4; 2.40797685e-9 8.48749924e-11])
end
