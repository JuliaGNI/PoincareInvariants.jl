@safetestset "getpaduanum, getdegree and checkpaduanum" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.Padua

    for n in [1:10..., 50, 123, 511, 10_000]
        paduanum = getpaduanum(n)
        @test paduanum == (n + 1) * (n + 2) ÷ 2
        @test getdegree(paduanum) == n

        @test_throws ArgumentError getdegree(paduanum + 1)
        @test_throws ArgumentError getdegree(paduanum - 1)

        @test nextpaduanum(paduanum - 1) == paduanum
        @test nextpaduanum(paduanum) == paduanum
        @test nextpaduanum(paduanum + 1) == getpaduanum(n + 1)
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

    @test Padua.fromcoeffsmat!(zeros(5, 4), reshape(1:20, 5, 4), 3) == [
        1  6  11 16;
        2  7  12  0;
        3  8   0  0;
        4  0   0  0;
        0  0   0  0
    ]
end

@safetestset "paduatransform! and invpaduatransform!" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.Padua
    using Statistics: mean
    using StaticArrays: SVector

    T0(x) = 1
	T1(x) = x
	T2(x) = 2x^2 - 1
	T3(x) = 4x^3 - 3x

    function T(v, cfs)
        cfs[1, 1] * T0(v[1]) * T0(v[2]) +

        cfs[1, 2] * T1(v[1]) * T0(v[2]) +
        cfs[2, 1] * T0(v[1]) * T1(v[2]) +

        cfs[1, 3] * T2(v[1]) * T0(v[2]) +
        cfs[2, 2] * T1(v[1]) * T1(v[2]) +
        cfs[3, 1] * T0(v[1]) * T2(v[2]) +

        cfs[1, 4] * T3(v[1]) * T0(v[2]) +
        cfs[2, 3] * T2(v[1]) * T1(v[2]) +
        cfs[3, 2] * T1(v[1]) * T2(v[2]) +
        cfs[4, 1] * T0(v[1]) * T3(v[2])
    end

    for degree in [4, 7, 21, 100]
        plan = PaduaTransformPlan{Float64}(degree)
        iplan = InvPaduaTransformPlan{Float64}(degree)
        points = getpaduapoints(degree)

        paduanum = getpaduanum(3)

        cfsmat = Padua.fromcoeffsmat!(zeros(3+2, 3+1), rand(3+2, 3+1), 3)
        cfsveclexfalse = Padua.fromcoeffsmat!(Vector{Float64}(undef, paduanum), cfsmat, 3, Val(false))
        cfsveclextrue = Padua.fromcoeffsmat!(Vector{Float64}(undef, paduanum), cfsmat, 3, Val(true))
        cfsarr = cat(cfsmat, cfsmat, cfsmat, cfsmat; dims=3)

        for i in 1:4
            cfsarr[1, 1, i] += i
        end

        vals = map(v -> T(v, cfsmat), points)
        valsmat = hcat(vals .+ 1, vals .+ 2, vals .+ 3, vals .+ 4)
        valsvecvec = SVector{4, Float64}.(eachrow(valsmat))

        vals2 = similar(vals)
        valsmat2 = similar(valsmat)

        cfsvec2 = Vector{Float64}(undef, getpaduanum(degree))
        cfsmat2 = zeros(degree+2, degree+1)
        cfsarr2 = zeros(degree+2, degree+1, 4)

        paduatransform!(cfsvec2, plan, vals, Val(true))
        @test mean(abs, cfsvec2[1:10] .- cfsveclextrue) < 5eps()
        @test mean(abs, cfsvec2[11:end]) < 5eps()

        paduatransform!(cfsvec2, plan, vals, Val(false))
        @test mean(abs, cfsvec2[1:10] .- cfsveclexfalse) < 5eps()
        @test mean(abs, cfsvec2[11:end]) < 5eps()

        paduatransform!(cfsmat2, plan, vals)
        @test mean(abs, cfsmat2[1:5, 1:4] .- cfsmat) < 5eps()
        @test mean(abs, cfsmat2[6:end, 1:end]) < 5eps()
        @test mean(abs, cfsmat2[1:end, 5:end]) < 5eps()

        invpaduatransform!(vals2, iplan, cfsmat2)
        @test mean(abs, vals2 .- vals) < 10eps()

        paduatransform!(cfsarr2, plan, valsmat)
        @test mean(abs, cfsarr2[1:5, 1:4, :] .- cfsarr) < 5eps()
        @test mean(abs, cfsarr2[6:end, 1:end, :]) < 5eps()
        @test mean(abs, cfsarr2[1:end, 5:end, :]) < 5eps()

        invpaduatransform!(valsmat2, iplan, cfsarr2)
        @test mean(abs, valsmat2 .- valsmat) < 10eps()

        paduatransform!(cfsarr2, plan, valsvecvec)
        @test mean(abs, cfsarr2[1:5, 1:4, :] .- cfsarr) < 5eps()
        @test mean(abs, cfsarr2[6:end, 1:end, :]) < 5eps()
        @test mean(abs, cfsarr2[1:end, 5:end, :]) < 5eps()
    end
end
