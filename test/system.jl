using VibronicToolkit: issimple, potential

using LinearAlgebra: Diagonal

let sys1 = DenseSystem(test_params...),
    sys2 = DenseSystem(test_params...),
    sys3 = DenseSystem(test_params_diag...),
    sys4 = DiagonalSystem(test_params_diag...),
    sys5 = DenseSystem(ones(1,1), [ones(1,1), ones(1,1,1), ones(1,1,1,1)]),
    sys6 = DenseSystem(ones(1,1), [ones(1,1), ones(1,1,1), ones(1,1,1,1), ones(1,1,1,1,1), ones(1,1,1,1,1,1)])

    @test sys1 == sys1
    @test sys1 == sys2
    @test sys1 != sys3
    @test sys1 != sys4
    @test sys1 != sys5
    @test sys1 != sys6

    @test sys2 == sys1
    @test sys2 == sys2
    @test sys2 != sys3
    @test sys2 != sys4
    @test sys2 != sys5
    @test sys2 != sys6

    @test sys3 != sys1
    @test sys3 != sys2
    @test sys3 == sys3
    @test sys3 == sys4
    @test sys3 != sys5
    @test sys3 != sys6

    @test sys4 != sys1
    @test sys4 != sys2
    @test sys4 == sys3
    @test sys4 == sys4
    @test sys4 != sys5
    @test sys4 != sys6

    @test sys5 != sys1
    @test sys5 != sys2
    @test sys5 != sys3
    @test sys5 != sys4
    @test sys5 == sys5
    @test sys5 != sys6

    @test sys6 != sys1
    @test sys6 != sys2
    @test sys6 != sys3
    @test sys6 != sys4
    @test sys6 != sys5
    @test sys6 == sys6
end

let sys1 = DenseSystem(test_freq, []),
    sys2 = DiagonalSystem(test_freq, [])

    @test sys1 == sys2
    @test issimple(sys1)
    @test sys1 == simplify(sys1)
    @test issimple(sys2)
    @test sys2 == simplify(sys2)
end

let sys1 = DenseSystem(test_params...),
    sys2 = DiagonalSystem(test_params_diag...),
    qs = [0.9, -0.23]

    expected_V = [-8.79136561e-2 -9.66581850e-1 -9.55687400e-1; -9.66581850e-1 3.65425699e-1 -9.44792950e-1; -9.55687400e-1 -9.44792950e-1 8.18765054e-1]
    @test isapprox(potential(sys1, qs), expected_V)
    @test isapprox(potential(sys2, qs), Diagonal(expected_V))
end

@testset "system_coef.jl" begin include("system_coef.jl") end
@testset "system_dense.jl" begin include("system_dense.jl") end
@testset "system_diagonal.jl" begin include("system_diagonal.jl") end
@testset "system_io.jl" begin include("system_io.jl") end
