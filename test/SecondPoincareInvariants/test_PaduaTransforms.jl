using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.PaduaTransforms

module ChebyshevTestUtils
    using LinearAlgebra: UpperTriangular

    export T0, T1, T2, T3, T6, T11
    export evalT, remove_lower_right

    T0(x) = evalpoly(x, (1,))
    T1(x) = evalpoly(x, (0, 1))
    T2(x) = evalpoly(x, (-1, 0, 2))
    T3(x) = evalpoly(x, (0, -3, 0, 4))
    T6(x) = evalpoly(x, (-1, 0, 18, 0, -48, 0, 32))
    T11(x) = evalpoly(x, (0, -11, 0, 220, 0, -1232, 0, 2816, 0, -2816, 0, 1024))

    function evalT(x::T, n::Integer) where T
        if n == 0
            return one(T)
        elseif n == 1
            return x
        elseif n > 1
            a = one(T)
            b = x
            for _ in 2:n
                b, a = 2x * b - a, b
            end
            return b
        else
            throw(ArgumentError("n must be greater than 0"))
        end
    end

    function evalT(x::T, y::T, coeffs::AbstractMatrix) where T
        axes(coeffs, 1) == axes(coeffs, 2) || error()
        s = zero(T)
        fi = firstindex(coeffs, 1)

        for i in axes(coeffs, 1)
            for j in axes(coeffs, 2)
                s += coeffs[i, j] * evalT(y, i - fi) * evalT(x, j - fi)
            end
        end

        s
    end

    remove_lower_right(A) = reverse(Matrix(UpperTriangular(A)); dims=2)
end  # ChebyshevTestUtils

@safetestset "ChebyshevTestUtils" begin
    using ..ChebyshevTestUtils

    @test evalT(0.225, 0) == 1.0
    @test evalT(0.225, 1) == 0.225
    @test evalT(0.225, 6) ≈ T6(0.225) atol=5eps()
    @test evalT(0.225, 11) == let x = 0.225
        1024x^11 - 2816x^9 + 2816x^7 - 1232x^5 + 220x^3 - 11x
    end
end

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

@safetestset "paduapoint and ispadua" begin
    using ..PaduaTransforms: paduapoint, ispadua

    @test paduapoint(Float64, 1, 3, 5)[1] ≈ cos(π * 1 / 5) atol=3eps()
    @test paduapoint(Float64, 1, 3, 5)[2] ≈ cos(π * 3 / 6) atol=3eps()

    @test paduapoint(Float64, 23, 42, 100)[1] ≈ cos(π * 23 / 100) atol=3eps()
    @test paduapoint(Float64, 23, 42, 100)[2] ≈ cos(π * 42 / 101) atol=3eps()

    @test eltype(paduapoint(Float32, 1, 2, 4)) == Float32
    @test eltype(paduapoint(Float64, 5, 6, 7)) == Float64

    @test ispadua(0, 0) == true
    @test ispadua(1, 0) == false
    @test ispadua(0, 1) == false
    @test ispadua(5, 4) == false
    @test ispadua(4, 4) == true
    @test ispadua(2, 0) == true
end

@safetestset "getpaduapoints" begin
    using ..PaduaTransforms

    @test getpaduapoints(1) isa NTuple{2, Vector{Float64}}
    @test getpaduapoints(Float32, 2) isa NTuple{2, Vector{Float32}}
    @test getpaduapoints(Float64, 3) isa NTuple{2, Vector{Float64}}

    @test getpaduapoints(1)[1] ≈ [1, 1, -1]
    @test getpaduapoints(1)[2] ≈ [1, -1, 0]

    @test getpaduapoints(Float32, 2)[1] ≈
        [cospi(0/2), cospi(0/2), cospi(1/2), cospi(1/2), cospi(2/2), cospi(2/2)]
    @test getpaduapoints(Float32, 2)[2] ≈
        [cospi(0/3), cospi(2/3), cospi(1/3), cospi(3/3), cospi(0/3), cospi(2/3)]


    @test getpaduapoints(Float32, 4)[1] ≈ [
        cospi(0/4), cospi(0/4), cospi(0/4),
        cospi(1/4), cospi(1/4), cospi(1/4),
        cospi(2/4), cospi(2/4), cospi(2/4),
        cospi(3/4), cospi(3/4), cospi(3/4),
        cospi(4/4), cospi(4/4), cospi(4/4),
    ]
    @test getpaduapoints(Float32, 4)[2] ≈ [
        cospi(0/5), cospi(2/5), cospi(4/5),
        cospi(1/5), cospi(3/5), cospi(5/5),
        cospi(0/5), cospi(2/5), cospi(4/5),
        cospi(1/5), cospi(3/5), cospi(5/5),
        cospi(0/5), cospi(2/5), cospi(4/5),
    ]

    @test getpaduapoints(Float32, 5)[1] ≈ [
        cospi(0/5), cospi(0/5), cospi(0/5), cospi(0/5),
        cospi(1/5), cospi(1/5), cospi(1/5),
        cospi(2/5), cospi(2/5), cospi(2/5), cospi(2/5),
        cospi(3/5), cospi(3/5), cospi(3/5),
        cospi(4/5), cospi(4/5), cospi(4/5), cospi(4/5),
        cospi(5/5), cospi(5/5), cospi(5/5),
    ]
    @test getpaduapoints(Float32, 5)[2] ≈ [
        cospi(0/6), cospi(2/6), cospi(4/6), cospi(6/6),
        cospi(1/6), cospi(3/6), cospi(5/6),
        cospi(0/6), cospi(2/6), cospi(4/6), cospi(6/6),
        cospi(1/6), cospi(3/6), cospi(5/6),
        cospi(0/6), cospi(2/6), cospi(4/6), cospi(6/6),
        cospi(1/6), cospi(3/6), cospi(5/6),
    ]

    @inferred getpaduapoints(3)
    @inferred getpaduapoints((x, y) -> x + y, 4)
    @inferred getpaduapoints((x, y) -> (x, y, x * x + y), 5)

    dopoints1 = getpaduapoints(11) do x, y
        x * y
    end

    @test dopoints1 isa Vector{Float64}
    @test length(dopoints1) == getpaduanum(11)
    @test dopoints1 ≈ getpaduapoints(11)[1] .* getpaduapoints(11)[2]

    dopoints4 = getpaduapoints(Float32, 5) do x, y
        x, y, x * y, x + y
    end

    @test dopoints4 isa NTuple{4, Vector{Float32}}
    @test all(dopoints4) do v
        length(v) == getpaduanum(5)
    end

    @test dopoints4[1] ≈ getpaduapoints(5)[1]
    @test dopoints4[2] ≈ getpaduapoints(5)[2]
    @test dopoints4[3] ≈ getpaduapoints(5)[1] .* getpaduapoints(5)[2]
    @test dopoints4[4] ≈ getpaduapoints(5)[1] .+ getpaduapoints(5)[2]
end

@safetestset "weight! and invweight!" begin
    using ..PaduaTransforms:
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
    using ..PaduaTransforms

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
    using ..PaduaTransforms

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

@safetestset "tocoeffsmat!" begin
    using ..PaduaTransforms

    mat = PaduaTransforms.tocoeffsmat!(zeros(5, 4), reshape(1:16, 4, 4))

    @test mat == [
        1.0  5.0   9.0  13.0
        2.0  6.0  10.0  14.0
        3.0  7.0  11.0  15.0
        4.0  8.0  12.0  16.0
        0.0  0.0   0.0   0.0
    ]
end

@safetestset "fromvalsmat!" begin
    using ..PaduaTransforms

    mat1 = reshape(1:20, 5, 4)
    out1 = zeros(10)

    PaduaTransforms.fromvalsmat!(out1, mat1, 3)

    @test out1 == [1, 3, 5,  7, 9,  11, 13, 15,  17, 19]

    mat2 = reshape(1:12, 4, 3)
    out2 = zeros(6)

    PaduaTransforms.fromvalsmat!(out2, mat2, 2)

    @test out2 == [1, 3,  6, 8,  9, 11]
end

@safetestset "1D paduatransform!" begin
    using ..PaduaTransforms
    using ..ChebyshevTestUtils

    @testset "Low degree tests" begin
        let n = 3; vals = getpaduapoints(n) do x, y
                1.5 * T0(x) * T1(y) + 2 * T3(x) * T0(y) + 100 * T2(x) * T1(y)
            end

            plan = PaduaTransformPlan{Float64}(n)

            out = paduatransform!(zeros(n+1, n+1), plan, vals)

            @test maximum(abs, out .- [
                  0 0   0 2
                1.5 0 100 0
                  0 0   0 0
                  0 0   0 0
            ]) / eps() < 100
        end

        let n = 6; vals = getpaduapoints(n) do x, y
                15 * T0(x) * T0(y) + 2.5 * T0(x) * T6(y) + 0.1 * T3(x) * T3(y)
            end

            plan = PaduaTransformPlan{Float64}(n)

            out = paduatransform!(zeros(n+1, n+1), plan, vals)

            @test maximum(abs, out .- [
                 15 0 0   0 0 0 0
                  0 0 0   0 0 0 0
                  0 0 0   0 0 0 0
                  0 0 0 0.1 0 0 0
                  0 0 0   0 0 0 0
                  0 0 0   0 0 0 0
                2.5 0 0   0 0 0 0
            ]) / eps() < 100
        end
    end

    @testset "Degree $n tests" for n in [11, 20, 29, 40, 51]
        plan = PaduaTransformPlan{Float64}(n)

        testcoeffs = remove_lower_right(PaduaTransforms.rand(n+1, n+1))

        vals = getpaduapoints(n) do x, y
            evalT(x, y, testcoeffs)
        end

        out = paduatransform!(zeros(n+1, n+1), plan, vals)

        @test maximum(abs, out .- testcoeffs) / eps() < 100
    end
end

@safetestset "1D invpaduatransform!" begin
    using ..PaduaTransforms
    using ..ChebyshevTestUtils

    @testset "Low degree tests" begin
        let n = 3; testvals = getpaduapoints(n) do x, y
                1.5 * T0(x) * T1(y) + 2 * T3(x) * T0(y) + 100 * T2(x) * T1(y)
            end

            invplan = InvPaduaTransformPlan{Float64}(n)

            coeffs = [
                  0 0   0 2;
                1.5 0 100 0;
                  0 0   0 0;
                  0 0   0 0
            ]

            vals = Vector{Float64}(undef, getpaduanum(n))
            invpaduatransform!(vals, invplan, coeffs)

            @test maximum(abs, vals .- testvals) / eps() < 100
        end

        let n = 6; testvals = getpaduapoints(n) do x, y
                15 * T0(x) * T0(y) + 2.5 * T0(x) * T6(y) + 0.1 * T3(x) * T3(y)
            end

            invplan = InvPaduaTransformPlan{Float64}(n)

            coeffs = [
                 15 0 0   0 0 0 0
                  0 0 0   0 0 0 0
                  0 0 0   0 0 0 0
                  0 0 0 0.1 0 0 0
                  0 0 0   0 0 0 0
                  0 0 0   0 0 0 0
                2.5 0 0   0 0 0 0
            ]

            vals = Vector{Float64}(undef, getpaduanum(n))
            invpaduatransform!(vals, invplan, coeffs)

            @test maximum(abs, vals .- testvals) / eps() < 100
        end
    end

    @testset "Degree $n tests" for n in [11, 20, 29, 40, 51]
        invplan = InvPaduaTransformPlan{Float64}(n)

        coeffs = rand(n+1, n+1)

        testvals = getpaduapoints(n) do x, y
            evalT(x, y, coeffs)
        end

        vals = Vector{Float64}(undef, getpaduanum(n))
        invpaduatransform!(vals, invplan, coeffs)

        @test maximum(abs, vals .- testvals) / eps() < 0.1 * n^3.3
    end
end

@safetestset "ND paduatransform! and invpaduatransform!" begin
    using ..PaduaTransforms

    n = 500
    D = 6
    plan = PaduaTransformPlan{Float64}(n)
    invplan = InvPaduaTransformPlan{Float64}(n)

    ndtestvals = rand(getpaduanum(n), D)
    ndcoeffs = zeros(n+1, n+1, D)
    paduatransform!(ndcoeffs, plan, ndtestvals)

    ndvals = Matrix{Float64}(undef, getpaduanum(n), D)
    invpaduatransform!(ndvals, invplan, ndcoeffs)

    @test maximum(abs, ndvals .- ndtestvals) / eps() < 100

    for i in 1:D
        coeffs = zeros(n+1, n+1)
        @test ndcoeffs[:, :, i] == paduatransform!(coeffs, plan, ndtestvals[:, i])
    end

    for i in 1:D
        vals = Vector{Float64}(undef, getpaduanum(n))
        @test ndvals[:, i] == invpaduatransform!(vals, invplan, ndcoeffs[:, :, i])
    end
end
