using VibronicToolkit: HamiltonianCoefficients, mode_indices

let cub = zeros(test_M, test_M, test_M, test_S, test_S),
    quart = ones(test_M, test_M, test_M, test_M, test_S, test_S),
    coef = HamiltonianCoefficients(test_params[2][1:3]..., cub, quart)

    @test keys(coef) == [0, 1, 2, 4]
    @test haskey(coef, 0)
    @test haskey(coef, 1)
    @test haskey(coef, 2)
    @test !haskey(coef, 3)
    @test haskey(coef, 4)
    @test !haskey(coef, 5)

    @test coef[0] == test_params[2][1]
    @test coef[1] == test_params[2][2]
    @test coef[2] == test_params[2][3]
    @test coef[3] == cub
    @test coef[4] == quart
    @test coef[5] == zeros(test_M, test_M, test_M, test_M, test_M, test_S, test_S)

    keys_collected = Int[]
    for (ord, val) in coef
        @test val == coef[ord]
        push!(keys_collected, ord)
    end
    @test keys_collected == keys(coef)

    @test length(mode_indices(coef[6])) == test_M^6

    @test (5.6 * coef)[4] == 5.6 * quart
end

let coef12 = HamiltonianCoefficients{1,2}(),
    coef21 = HamiltonianCoefficients{2,1}()

    @test coef12[0] == zeros(1, 1)
    @test coef12[1] == zeros(2, 1, 1)
    @test coef12[2] == zeros(2, 2, 1, 1)
    @test coef12[3] == zeros(2, 2, 2, 1, 1)

    @test coef21[0] == zeros(2, 2)
    @test coef21[1] == zeros(1, 2, 2)
    @test coef21[2] == zeros(1, 1, 2, 2)
    @test coef21[3] == zeros(1, 1, 1, 2, 2)
end

# Not enough information to determine the number of modes.
@test_throws DomainError HamiltonianCoefficients(ones(3, 3))
# Surface mismatch.
@test_throws DomainError HamiltonianCoefficients{2,1}(ones(1, 1))
@test_throws DomainError HamiltonianCoefficients{2,1}(ones(1, 1), ones(1, 1, 1))
@test_throws DomainError HamiltonianCoefficients(ones(1, 1), ones(1, 2, 2))
# Mode mismatch.
@test_throws DomainError HamiltonianCoefficients{1,2}(ones(1, 1), ones(1, 1, 1))
@test_throws DomainError HamiltonianCoefficients(ones(1, 1), ones(1, 1, 1), ones(2, 2, 1, 1))
