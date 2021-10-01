@safetestset "getpaduanum, getdegree and checkpaduanum" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.Padua

    for n in [1:10..., 50, 123, 511, 10_000]
        paduanum = getpaduanum(n)
        @test paduanum == (n + 1) * (n + 2) ÷ 2
        @test isinteger(getdegree(paduanum))
        @test (checkpaduanum(paduanum); true)

        @test_throws ArgumentError checkpaduanum(paduanum + 1)
        @test_throws ArgumentError checkpaduanum(paduanum - 1)
    end
end

@safetestset "chebyshevpoint and paduapoint" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.Padua:
        chebyshevpoint, paduapoint
    
    @test chebyshevpoint(Float64, 0, 0, 1, 1) ≈ [ 1,  1]
    @test chebyshevpoint(Float64, 0, 1, 1, 1) ≈ [ 1, -1]
    @test chebyshevpoint(Float64, 1, 0, 1, 1) ≈ [-1,  1]
    @test chebyshevpoint(Float64, 1, 1, 1, 1) ≈ [-1, -1]
    
    @test chebyshevpoint(Float64, 1, 3, 2, 5) ≈ [cos(π * 1 / 2), cos(π * 3 / 5)] atol=3eps()

    @test eltype(chebyshevpoint(Float32, 1, 2, 3, 4)) == Float32
    @test eltype(chebyshevpoint(Float64, 5, 6, 7, 8)) == Float64

    @test paduapoint(Float64, 1, 2, 3) === chebyshevpoint(Float64, 1, 2, 3, 3+1)
    @test paduapoint(Float32, 12, 73, 150) === chebyshevpoint(Float32, 12, 73, 150, 150+1)

    @test [paduapoint(Float32, x, y, 1) for y in 0:1+1, x in 0:0] ≈ [
        [1,  1],
        [1,  0],
        [1, -1]
    ]

    @test [paduapoint(Float64, x, y, 1) for y in 0:1+1, x in 1:1] ≈ [
        [-1,  1],
        [-1,  0],
        [-1, -1]
    ]
end

@safetestset "ispadua and getpaduapoints" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.Padua

    @test Padua.ispadua(0, 0) == true
    @test Padua.ispadua(1, 0) == false
    @test Padua.ispadua(0, 1) == false
    @test Padua.ispadua(5, 4) == false
    @test Padua.ispadua(4, 4) == true
    @test Padua.ispadua(2, 0) == true

    @test getpaduapoints(1) ≈ [
        [1, 1],
        [1, -1],
        [-1, 0]
    ]

    for n in 1:20
        pnts = getpaduapoints(n)
        @test pnts == [Padua.paduapoint(Float64, x, y, n) for y in 0:n+1, x in 0:n if Padua.ispadua(x, y)]

        pnts32 = getpaduapoints(Float32, n)
        @test pnts32 == [Padua.paduapoint(Float32, x, y, n) for y in 0:n+1, x in 0:n if Padua.ispadua(x, y)]
    end
end

@safetestset "weight!" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.Padua: weight!

    @test weight!(ones(4+2, 4+1), 4) == [
        0.025  0.05  0.05  0.05  0.025;
         0.05   0.1   0.1   0.1   0.05;
         0.05   0.1   0.1   0.1   0.05;
         0.05   0.1   0.1   0.1   0.05;
         0.05   0.1   0.1   0.1   0.05;
        0.025  0.05  0.05  0.05  0.025
    ]
end

@safetestset "tovalsmat!" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.Padua

    @test Padua.tovalsmat!(ones(3 + 2, 3 + 1), 1:getpaduanum(3), 3) == [
        1 0 6 0 ;
        0 4 0 9 ;
        2 0 7 0 ;
        0 5 0 10;
        3 0 8 0
    ]

    @test Padua.tovalsmat!(ones(2 + 2, 2 + 1), 1:getpaduanum(2), 2) == [
        1 0 5;
        0 3 0;
        2 0 6;
        0 4 0;
    ]
end

@safetestset "fromcoeffsmat!" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.Padua

    mat = [(x, y) for y in 0:2+1, x in 0:2]
    to = similar(mat, getpaduanum(2))

    Padua.fromcoeffsmat!(to, mat, 2, Val(true))
    @test to == [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]

    Padua.fromcoeffsmat!(to, mat, 2, Val(false))
    @test to == [(0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0)]
end

@safetestset "paduatransform! and PaduaTransformPlan" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.Padua

    T0(x) = 1
	T1(x) = x
	T2(x) = 2x^2 - 1
	T3(x) = 4x^3 - 3x

    function T(v, cfs)
        cfs[1] * T0(v[1]) * T0(v[2]) +

        cfs[2] * T1(v[1]) * T0(v[2]) +
        cfs[3] * T0(v[1]) * T1(v[2]) +

        cfs[4] * T2(v[1]) * T0(v[2]) +
        cfs[5] * T1(v[1]) * T1(v[2]) +
        cfs[6] * T0(v[1]) * T2(v[2]) +

        cfs[7] * T3(v[1]) * T0(v[2]) +
        cfs[8] * T2(v[1]) * T1(v[2]) +
        cfs[9] * T1(v[1]) * T2(v[2]) +
        cfs[10] * T0(v[1]) * T3(v[2])
    end

    for degree in 4:10
        plan = PaduaTransformPlan{Float64}(degree)
        points = getpaduapoints(degree)

        coeffs = rand(10)
        vals = map(v -> T(v, coeffs), points)

        out = Vector{Float64}(undef, getpaduanum(degree))

        paduatransform!(out, plan, vals, Val(true))

        @test out[1:10] ≈ coeffs atol=10eps()
        @test out[11:end] ≈ zeros(length(out) - 10) atol=10eps()
    end
end

@safetestset "multi-dimensional paduatransform!" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.Padua
    using StaticArrays

    degree = 20
    dims = 5
    plan = PaduaTransformPlan{Float64}(degree)

    valsmat = rand(getpaduanum(degree), dims)
    outmat = Matrix{Float64}(undef, getpaduanum(degree), dims)

    paduatransform!(outmat, plan, valsmat, Val(false))

    valsvec = map(SVector{dims, Float64}, eachrow(valsmat))
    outvec = Matrix{Float64}(undef, getpaduanum(degree), dims)

    paduatransform!(outvec, plan, valsvec, Val(false))

    outcompare = Matrix{Float64}(undef, getpaduanum(degree), dims)

    for i in 1:dims
        paduatransform!(view(outcompare, :, i), plan, valsmat[:, i], Val(false))
    end

    @test outcompare == outmat
    @test outcompare == outvec
end
