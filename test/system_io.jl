let sys = read("data/dense.json", DenseSystem)

    @test isapprox(sys, DenseSystem(test_params_flat...))

    buf = IOBuffer()
    write(buf, sys)
    seekstart(buf)
    @test read(buf, DenseSystem) == sys
end

let sys = read("data/diagonal.json", DenseSystem)

    @test isapprox(sys, DenseSystem(test_params_flat_diag...))

    buf = IOBuffer()
    write(buf, sys)
    seekstart(buf)
    @test read(buf, DenseSystem) == sys
end

let sys = read("data/diagonal.json", DiagonalSystem)

    @test isapprox(sys, DiagonalSystem(test_params_flat_diag...))

    buf = IOBuffer()
    write(buf, sys)
    seekstart(buf)
    @test read(buf, DiagonalSystem) == sys
end

@test_throws SurfaceCouplingException read("data/dense.json", DiagonalSystem)
