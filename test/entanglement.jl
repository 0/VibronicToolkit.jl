using VibronicToolkit: S_renyi, S_vn

using LinearAlgebra: tr

let rho = [0.75 0.2; 0.2 0.25]

    @test S_vn(rho) == -tr(rho * log(rho))
    @test S_renyi(rho) == -log(tr(rho^2))
end
