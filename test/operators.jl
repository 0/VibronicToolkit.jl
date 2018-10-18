using VibronicToolkit: a, Basis, id, mkid, mkop, n, operators, q,
                       trial_uniform, vectors

using LinearAlgebra: I

let sys = DenseSystem(test_params_trivial...),
    size = 3,
    basis = Basis(sys, size)

    @test id(basis) == I
    @test isapprox(a(basis)' * a(basis), n(basis))

    q_ = q(basis)
    ip_ = -(a(basis)' - a(basis)) / sqrt(2.0)
    ho = (-ip_^2 + q_^2)/2
    @test isapprox(ho[1:(end-1),1:(end-1)], n(basis)[1:(end-1),1:(end-1)]+0.5I)

    @test mkid(basis) == I
    @test mkop(basis, a, 1) == kron(a(basis), id(basis))
    @test mkop(basis, q, 2) == kron(id(basis), q(basis))

    h0, V = operators(basis, sys)
    h0_exp = zeros(basis.dim, basis.dim)
    for s in 1:test_S
        idx = ((s-1)*basis.dim1+1):(s*basis.dim1)
        h0_exp[idx, idx] .+= I + kron(n(basis)+0.5I, id(basis)) + kron(id(basis), n(basis)+0.5I)
    end
    @test h0 == h0_exp
    V_exp = zeros(basis.dim, basis.dim)
    for s1 in 1:test_S
        for s2 in 1:test_S
            s1 == s2 && continue
            idx1 = ((s1-1)*basis.dim1+1):(s1*basis.dim1)
            idx2 = ((s2-1)*basis.dim1+1):(s2*basis.dim1)
            V_exp[idx2, idx1] .+= Matrix(I, basis.dim1, basis.dim1)
        end
    end
    @test V == V_exp

    Vs = vectors(basis)
    conf1 = Vs[:, 24]
    @test conf1 == [1, 2, 3]
    vector1 = zeros(basis.dim)
    vector1[24] = 1.0
    conf2 = Vs[:, 25]
    @test conf2 == [2, 0, 3]
    vector2 = zeros(basis.dim)
    vector2[25] = sqrt(2.0)^2
    @test vector2 == kron(Matrix(I, test_S, test_S),
                          mkop(basis, a, 2)^2 * mkop(basis, a, 1)') * vector1

    @test isapprox(trial_uniform(basis), [3.54490770, 0.0, 2.50662827, 0.0,
                                          0.0, 0.0, 2.50662827, 0.0,
                                          1.77245385,
                                          3.54490770, 0.0, 2.50662827, 0.0,
                                          0.0, 0.0, 2.50662827, 0.0,
                                          1.77245385,
                                          3.54490770, 0.0, 2.50662827, 0.0,
                                          0.0, 0.0, 2.50662827, 0.0,
                                          1.77245385])
end
