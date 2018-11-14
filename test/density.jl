let sys = DenseSystem(test_params...),
    beta = 12.34,
    basis_size = 10,
    density = diagonal_density(sys, beta, basis_size, (-1.0, 2.0, -3.0, 4.0); lengths=(2, 3), progress_output=devnull)

    @test isapprox(density[:, :, 1, 1], [8.430118687e-6 4.249521714e-7; 5.397109667e-2 2.574071907e-3; 7.428088944e-9 3.353417148e-10])
    @test isapprox(density[:, :, 2, 3], [2.829791299e-6 1.434083980e-7; 1.818544383e-2 8.602008088e-4; 2.466985602e-9 1.088707521e-10])
end

let sys = DenseSystem(test_params...),
    trial = UniformTrialWavefunction(sys),
    beta = 12.34,
    basis_size = 10,
    (density, density_exact) = diagonal_density_pigs(sys, trial, beta, basis_size, (-1.0, 2.0, -3.0, 4.0); lengths=(2, 3), progress_output=devnull)

    @test isapprox(density[:, :, 1, 1], [8.430118481e-6 4.249521587e-7; 5.397109535e-2 2.574071826e-3; 7.428088751e-9 3.353417032e-10])
    @test isapprox(density[:, :, 2, 3], [2.829791400e-6 1.434084025e-7; 1.818544447e-2 8.602008361e-4; 2.466985689e-9 1.088707557e-10])
    @test isapprox(density_exact[:, :, 1, 1], [8.430118687e-6 4.249521714e-7; 5.397109667e-2 2.574071907e-3; 7.428088944e-9 3.353417148e-10])
    @test isapprox(density_exact[:, :, 2, 3], [2.829791299e-6 1.434083980e-7; 1.818544383e-2 8.602008088e-4; 2.466985602e-9 1.088707521e-10])
end
