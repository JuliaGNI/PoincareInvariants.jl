using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.PaduaTransforms

@safetestset "getpaduanum, getdegree and nextpaduanum" begin
    using ..PaduaTransforms

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

@safetestset "paduapoint, ispadua and getpaduapoints" begin
    using ..PaduaTransforms: ispadua, paduapoint
    using ..PaduaTransforms

    @test paduapoint(Float64, 1, 3, 5)[1] ≈ cos(π * 1 / 5) atol=3eps()
    @test paduapoint(Float64, 1, 3, 5)[2] ≈ cos(π * 3 / 6) atol=3eps()

    @test paduapoint(Float64, 23, 42, 100)[1] ≈ cos(π * 23 / 100) atol=3eps()
    @test paduapoint(Float64, 23, 42, 100)[2] ≈ cos(π * 42 / 101) atol=3eps()

    @test eltype(paduapoint(Float32, 1, 2, 4)) == Float32
    @test eltype(paduapoint(Float64, 5, 6, 7)) == Float64

    # collect tuple to vector, so ≈ works
    @test [collect(paduapoint(Float32, x, y, 1)) for y in 0:1+1, x in 0:0] ≈ [
        [1,  1],
        [1,  0],
        [1, -1]
    ]

    @test ispadua(0, 0) == true
    @test ispadua(1, 0) == false
    @test ispadua(0, 1) == false
    @test ispadua(5, 4) == false
    @test ispadua(4, 4) == true
    @test ispadua(2, 0) == true

    @test getpaduapoints(1) ≈ [
         1  1;
         1 -1;
        -1  0
    ]

    @test getpaduapoints(4) ≈ [
        cospi(0/4) cospi(0/5);
        cospi(0/4) cospi(2/5);
        cospi(0/4) cospi(4/5);
        cospi(1/4) cospi(1/5);
        cospi(1/4) cospi(3/5);
        cospi(1/4) cospi(5/5);
        cospi(2/4) cospi(0/5);
        cospi(2/4) cospi(2/5);
        cospi(2/4) cospi(4/5);
        cospi(3/4) cospi(1/5);
        cospi(3/4) cospi(3/5);
        cospi(3/4) cospi(5/5);
        cospi(4/4) cospi(0/5);
        cospi(4/4) cospi(2/5);
        cospi(4/4) cospi(4/5)
    ]

    torowvec(t) = [t[1] t[2]]
    for n in [21, 22, 23, 100]
        pnts = getpaduapoints(n)

        @test size(pnts) == (getpaduanum(n), 2)

        testpnts = reduce(vcat, [
            paduapoint(Float64, x, y, n) |> torowvec
            for y in 0:n+1, x in 0:n if ispadua(x, y)
        ])

        @test pnts == testpnts
    end

    dopoints = getpaduapoints(5) do x, y
        x, y, x * y, x + y
    end

    @test size(dopoints) == (getpaduanum(5), 4)

    @test dopoints[:, 1] ≈ getpaduapoints(5)[:, 1]
    @test dopoints[:, 2] ≈ getpaduapoints(5)[:, 2]
    @test dopoints[:, 3] ≈ getpaduapoints(5)[:, 1] .* getpaduapoints(5)[:, 2]
    @test dopoints[:, 4] ≈ getpaduapoints(5)[:, 1] .+ getpaduapoints(5)[:, 2]
end

@safetestset "weight! and invweight!" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.PaduaTransforms:
        weight!, invweight!

    @test weight!(ones(4+2, 4+1), 4) == [
        0.025  0.05  0.05  0.05  0.025;
         0.05   0.1   0.1   0.1   0.05;
         0.05   0.1   0.1   0.1   0.05;
         0.05   0.1   0.1   0.1   0.05;
         0.05   0.1   0.1   0.1   0.05;
        0.025  0.05  0.05  0.05  0.025
    ]

    @test invweight!(ones(4+2, 4+1)) == [
          1    0.5    0.5    0.5    1;
        0.5   0.25   0.25   0.25  0.5;
        0.5   0.25   0.25   0.25  0.5;
        0.5   0.25   0.25   0.25  0.5;
        0.5   0.25   0.25   0.25  0.5;
          1    0.5    0.5    0.5    1
    ]
end

@safetestset "tovalsmat!" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.PaduaTransforms

    @test PaduaTransforms.tovalsmat!(ones(3 + 2, 3 + 1), 1:getpaduanum(3), 3) == [
        1 0 6 0 ;
        0 4 0 9 ;
        2 0 7 0 ;
        0 5 0 10;
        3 0 8 0
    ]

    @test PaduaTransforms.tovalsmat!(ones(2 + 2, 2 + 1), 1:getpaduanum(2), 2) == [
        1 0 5;
        0 3 0;
        2 0 6;
        0 4 0;
    ]
end

@safetestset "fromcoeffsmat!" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.PaduaTransforms

    mat = [(x, y) for y in 0:2+1, x in 0:2]
    to = similar(mat, getpaduanum(2))

    PaduaTransforms.fromcoeffsmat!(to, mat, 2, Val(true))
    @test to == [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]

    PaduaTransforms.fromcoeffsmat!(to, mat, 2, Val(false))
    @test to == [(0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0)]

    @test PaduaTransforms.fromcoeffsmat!(zeros(4, 4), reshape(1:20, 5, 4), 3) == [
        1  6  11 16;
        2  7  12  0;
        3  8   0  0;
        4  0   0  0
    ]
end

@safetestset "paduatransform! and invpaduatransform!" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.PaduaTransforms
    using Statistics: mean
    using StaticArrays: SVector

    T0(x) = 1
	T1(x) = x
	T2(x) = 2x^2 - 1
	T3(x) = 4x^3 - 3x

    for degree in [4, 7, 21, 100]
        plan = PaduaTransformPlan{Float64}(degree)
        iplan = InvPaduaTransformPlan{Float64}(degree)
        points = getpaduapoints(degree)

        paduanum = getpaduanum(3)

        randcfs = rand(3+2, 3+1)

        cfsmat = PaduaTransforms.fromcoeffsmat!(zeros(3+1, 3+1), randcfs, 3)
        cfsveclexfalse = PaduaTransforms.fromcoeffsmat!(Vector{Float64}(undef, paduanum), randcfs, 3, Val(false))
        cfsveclextrue = PaduaTransforms.fromcoeffsmat!(Vector{Float64}(undef, paduanum), randcfs, 3, Val(true))
        cfsarr = cat(cfsmat, cfsmat, cfsmat, cfsmat; dims=3)

        for i in 1:4
            cfsarr[1, 1, i] += i
        end

        vals = getpaduapoints(degree) do x, y
            cfsmat[1, 1] * T0(x) * T0(y) +

            cfsmat[1, 2] * T1(x) * T0(y) +
            cfsmat[2, 1] * T0(x) * T1(y) +

            cfsmat[1, 3] * T2(x) * T0(y) +
            cfsmat[2, 2] * T1(x) * T1(y) +
            cfsmat[3, 1] * T0(x) * T2(y) +

            cfsmat[1, 4] * T3(x) * T0(y) +
            cfsmat[2, 3] * T2(x) * T1(y) +
            cfsmat[3, 2] * T1(x) * T2(y) +
            cfsmat[4, 1] * T0(x) * T3(y)
        end

        valsmat = hcat(vals .+ 1, vals .+ 2, vals .+ 3, vals .+ 4)
        valsvecvec = SVector{4, Float64}.(eachrow(valsmat))

        vals2 = similar(vals)
        valsmat2 = similar(valsmat)

        cfsvec2 = Vector{Float64}(undef, getpaduanum(degree))
        cfsmat2 = zeros(degree+1, degree+1)
        cfsarr2 = zeros(degree+1, degree+1, 4)

        paduatransform!(cfsvec2, plan, vals, Val(true))
        @test mean(abs, cfsvec2[1:10] .- cfsveclextrue) < 5eps()
        @test mean(abs, cfsvec2[11:end]) < 5eps()

        paduatransform!(cfsvec2, plan, vals, Val(false))
        @test mean(abs, cfsvec2[1:10] .- cfsveclexfalse) < 5eps()
        @test mean(abs, cfsvec2[11:end]) < 5eps()

        paduatransform!(cfsmat2, plan, vals)
        @test mean(abs, cfsmat2[1:4, 1:4] .- cfsmat) < 5eps()
        @test mean(abs, cfsmat2[5:end, 1:end]) < 5eps()
        @test mean(abs, cfsmat2[1:end, 5:end]) < 5eps()

        invpaduatransform!(vals2, iplan, cfsmat2)
        @test mean(abs, vals2 .- vals) < 10eps()

        paduatransform!(cfsarr2, plan, valsmat)
        @test mean(abs, cfsarr2[1:4, 1:4, :] .- cfsarr) < 5eps()
        @test mean(abs, cfsarr2[5:end, 1:end, :]) < 5eps()
        @test mean(abs, cfsarr2[1:end, 5:end, :]) < 5eps()

        invpaduatransform!(valsmat2, iplan, cfsarr2)
        @test mean(abs, valsmat2 .- valsmat) < 10eps()

        paduatransform!(cfsarr2, plan, valsvecvec)
        @test mean(abs, cfsarr2[1:4, 1:4, :] .- cfsarr) < 5eps()
        @test mean(abs, cfsarr2[5:end, 1:end, :]) < 5eps()
        @test mean(abs, cfsarr2[1:end, 5:end, :]) < 5eps()
    end
end
