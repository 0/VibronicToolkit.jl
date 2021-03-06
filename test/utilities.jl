using VibronicToolkit: @cat

let A = collect(reshape(1:9, 3, 3))

    @test diag(A, 1, 2) == diag(A)
end

let A = collect(reshape(1:8, 2, 2, 2))

    @test diag(A, 1, 2) == [1 5; 4 8]
    @test diag(A, 1, 3) == [1 3; 6 8]
    @test diag(A, 2, 3) == [1 7; 2 8]
end

let A = collect(reshape(1:72, 2, 3, 2, 3, 2)),
    B = diag(A, 2, 4)

    for i in 1:2, j in 1:3, k in 1:2, l in 1:2
        @test B[i, j, k, l] == A[i, j, k, j, l]
    end
    C = diag(diag(A, 1, 3), 1, 4)
    for i in 1:2, j in 1:3, k in 1:3
        @test C[i, j, k] == A[i, j, i, k, i]
    end
    @test_throws DomainError diag(A, 1, 2)
end

let A = collect(reshape(1:72, 2, 3, 2, 3, 2)),
    B = collect(reshape(73:180, 2, 3, 2, 3, 3)),
    C = collect(reshape(181:216, 2, 3, 2, 3))

    @test A == reshape(1:72, 2, 3, 2, 3, 2)
    @cat(A, B)
    @test A == reshape(1:180, 2, 3, 2, 3, 5)
    @cat(A, C)
    @test A == reshape(1:216, 2, 3, 2, 3, 6)
end
