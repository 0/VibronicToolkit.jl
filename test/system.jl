include("system_dense.jl")
include("system_diagonal.jl")
include("system_io.jl")

let sys1 = DenseSystem(test_params...),
    sys2 = DenseSystem(test_params...),
    sys3 = DenseSystem(test_params_diag...),
    sys4 = DiagonalSystem(test_params_diag...),
    sys5 = DenseSystem(ones(1,1), ones(1,1), ones(1,1,1), ones(1,1,1,1))

    @test sys1 == sys1
    @test sys1 == sys2
    @test sys1 != sys3
    @test sys1 != sys4
    @test sys1 != sys5

    @test sys2 == sys1
    @test sys2 == sys2
    @test sys2 != sys3
    @test sys2 != sys4
    @test sys2 != sys5

    @test sys3 != sys1
    @test sys3 != sys2
    @test sys3 == sys3
    @test sys3 == sys4
    @test sys3 != sys5

    @test sys4 != sys1
    @test sys4 != sys2
    @test sys4 == sys3
    @test sys4 == sys4
    @test sys4 != sys5

    @test sys5 != sys1
    @test sys5 != sys2
    @test sys5 != sys3
    @test sys5 != sys4
    @test sys5 == sys5
end
