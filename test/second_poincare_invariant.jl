using FastTransforms: PaduaTransformPlan, paduapoints

@testset "Coefficient and Point Counts" begin
    for n in 1:11:10_000
        padua_num = PoincareInvariants.get_padua_num(n)
        @test padua_num == (n + 1) * (n + 2) ÷ 2
        @test isinteger(PoincareInvariants.get_n(padua_num))
        @test PoincareInvariants.get_n(padua_num) == n
        @test (PoincareInvariants.check_padua_num(padua_num); true)
        @test_throws ArgumentError PoincareInvariants.check_padua_num(padua_num + 1)
        @test_throws ArgumentError PoincareInvariants.check_padua_num(padua_num - 1)
        
        @test PoincareInvariants.get_uu_coeff_num(n) == n * (n + 1) ÷ 2
        @test PoincareInvariants.get_uu_point_num(n) == n^2
        @test PoincareInvariants.get_uu_coeff_num(n) == (PoincareInvariants.get_uu_point_num(n) - n) ÷ 2 + n
    end
end

@testset "PoincareInvariant2 Constructor" begin
    @test PoincareInvariant2 <: AbstractPoincareInvariant

    N = 2
    T = Float64

    padua_num = 10_011
    n = PoincareInvariants.get_n(padua_num)
    ppoints = get_padua_points(padua_num)
    pinv = PoincareInvariant2{N, T}(padua_num::Integer)

    @test pinv.n == n

    @test pinv.cc_coeffs isa Matrix{Float64}
    @test size(pinv.cc_coeffs) == (padua_num, N)

    @test size(pinv.D1toUU) == (PoincareInvariants.get_uu_coeff_num(n), padua_num)
    @test size(pinv.D2toUU) == (PoincareInvariants.get_uu_coeff_num(n), padua_num)
    @test size(pinv.CCtoUU) == (PoincareInvariants.get_uu_coeff_num(n), padua_num)
end

@testset "1D Free Particle" begin
    ppoints = paduapoints(140)

    pinv = PoincareInvariant2{2, Float64}(10_011)

    function free_particle!(init, δt)
        init[1] += init[2] .* δt
    end

    Ω = [0 -1; 1 0]

    @test compute(pinv, ppoints, Ω) ≈ 4

    free_particle!.(eachrow(ppoints), 10)
    @test compute(pinv, ppoints, Ω) ≈ 4

    free_particle!.(eachrow(ppoints), 10)
    @test compute(pinv, ppoints, Ω) ≈ 4
end
