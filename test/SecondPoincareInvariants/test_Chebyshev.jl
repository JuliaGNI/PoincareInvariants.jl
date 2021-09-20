@safetestset "Coefficient and Point Counts" begin
    using PoincareInvariants.SecondPoincareInvariants: getpaduanum, checkpaduanum,
        nextpaduanum, getdegree, getcoeffnum, getpointnum

    for n in [1:100..., 135, 752, 1000, 5531]
        paduanum = getpaduanum(n)
        @test paduanum == (n + 1) * (n + 2) ÷ 2
        @test isinteger(getdegree(paduanum))
        @test getdegree(paduanum) == n
        @test (checkpaduanum(paduanum); true)
        @test_throws ArgumentError checkpaduanum(paduanum + 1)
        @test_throws ArgumentError checkpaduanum(paduanum - 1)

        @test nextpaduanum(paduanum) == paduanum
        @test nextpaduanum(paduanum - 1) == paduanum
        @test nextpaduanum(paduanum + 1) == getpaduanum(n + 1)
        
        @test getcoeffnum(n) == paduanum
        @test getcoeffnum(n - 1) == n * (n + 1) ÷ 2
        @test getpointnum(n - 1) == n^2
        @test getcoeffnum(n - 1) == (getpointnum(n - 1) - n) ÷ 2 + n
    end
end

@safetestset "getpaduapoints" begin
    using PoincareInvariants.SecondPoincareInvariants: getpaduapoints, getpaduanum
    using FastTransforms: paduapoints

    @test eltype(getpaduapoints(21)) == Float64

    @testset "getpaduapoints($T, $n)" for T in [Float32, Float64], n in [1:25..., 50, 75, 100]
        N = getpaduanum(n)
        ftpoints = (paduapoints(T, n) .+ 1) ./ 2
        pipoints = getpaduapoints(T, n)
        
        @test ftpoints ≈ pipoints
        @test eltype(pipoints) == T
    end
end

@safetestset "paduatransform!" begin
    using PoincareInvariants.SecondPoincareInvariants: getpaduapoints, paduatransform!
    using FastTransforms: plan_paduatransform!

    f(v) = v[1] * sin(v[2])

    paduapoints = getpaduapoints(20)
    v = f.(eachrow(paduapoints))

    plan = plan_paduatransform!(v, Val{false})

    coeffs = plan * copy(v)

    out = similar(v)
    paduatransform!(out, v, plan)

    @test out == coeffs
    @test v == f.(eachrow(paduapoints))
end
