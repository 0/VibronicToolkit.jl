using VibronicToolkit: a, Basis, id, mkid, mkop, n, operators, q

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
end
