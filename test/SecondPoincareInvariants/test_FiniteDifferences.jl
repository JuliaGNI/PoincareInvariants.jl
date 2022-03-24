using PoincareInvariants.SecondPoincareInvariants.FiniteDifferences

@safetestset "getpoints, getpointspec and getpointnum" begin
    using ..FiniteDifferences: FiniteDiffPlan, getpoints, getpointnum, getpointspec

    @test getpointspec(50, FiniteDiffPlan) == (9, 9)
    @test getpointspec(87, FiniteDiffPlan) == (11, 11)
    @test getpointspec((53, 42), FiniteDiffPlan) == (53, 43)

    @test getpointnum((3, 5)) == 15
    @test getpointnum((423789, 326)) == 423789 * 326
    @test getpointnum((123, 908343)) == 123 * 908343

    @inferred NTuple{3, Vector{Float64}} getpoints((x, y) -> (x, y, x+y), Float64, (7, 5), FiniteDiffPlan)
    @inferred Vector{Float32} getpoints((x, y) -> sin(x), Float32, (9, 9), FiniteDiffPlan)

    @test getpoints((x, y) -> (x, y), Float64, (5, 3), FiniteDiffPlan)[1] ≈ [
        0   , 0   , 0   ,
        0.25, 0.25, 0.25,
        0.5 , 0.5 , 0.5 ,
        0.75, 0.75, 0.75,
        1.0 , 1.0 , 1.0 ]
    @test getpoints((x, y) -> (x, y), Float64, (5, 3), FiniteDiffPlan)[2] ≈ [
        0  , 0.5, 1.0,
        0  , 0.5, 1.0,
        0  , 0.5, 1.0,
        0  , 0.5, 1.0,
        0  , 0.5, 1.0]

    reshapepnts = getpoints((x, y) -> (4x, 6y), Float64, (5, 7), FiniteDiffPlan)
    @test reshape(reshapepnts[1], 7, 5) ≈ [
        0 1 2 3 4;
        0 1 2 3 4;
        0 1 2 3 4;
        0 1 2 3 4;
        0 1 2 3 4;
        0 1 2 3 4;
        0 1 2 3 4]
    @test reshape(reshapepnts[2], 7, 5) ≈ [
        0 0 0 0 0;
        1 1 1 1 1;
        2 2 2 2 2;
        3 3 3 3 3;
        4 4 4 4 4;
        5 5 5 5 5;
        6 6 6 6 6]

    @test getpoints((x, y) -> y, Float64, (11, 11), FiniteDiffPlan) ≈ Float64[
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    ]

    ps = getpointspec(6342, FiniteDiffPlan)
    testpnts = getpoints((x, y) -> (x, y), Float64, ps, FiniteDiffPlan)
    pnts4 = getpoints(Float64, ps, FiniteDiffPlan) do x, y
        x, y, x * y, x + y
    end
    @test pnts4 isa NTuple{4, Vector{Float64}}
    @test all(pnts4) do v
        length(v) == getpointnum(ps)
    end

    @test pnts4[1] ≈ testpnts[1]
    @test pnts4[2] ≈ testpnts[2]
    @test pnts4[3] ≈ testpnts[1] .* testpnts[2]
    @test pnts4[4] ≈ testpnts[1] .+ testpnts[2]
end

@safetestset "Differentiation" begin
    using ..FiniteDifferences: getpoints, FiniteDiffPlan, differentiate!

    # Finite difference derivatives should be exact for a quadratic
    let f(x, y)  = 1 + 0.75*x + y - 0.5*x*x - y*y - 1.5*x*y
        fx(x, y) = 0.75 - x - 1.5*y
        fy(x, y) = 1 - 2*y - 1.5*x

        # small
        sdims = (3, 3)
        svals = getpoints(f, Float64, sdims, FiniteDiffPlan)

        s∂x = zeros(Float64, length(svals))
        s∂y = zeros(Float64, length(svals))

        differentiate!(s∂x, s∂y, svals, sdims)

        sfx = getpoints(fx, Float64, sdims, FiniteDiffPlan)
        @test maximum(abs, s∂x .- sfx) / eps() < 50

        sfy = getpoints(fy, Float64, sdims, FiniteDiffPlan)
        @test maximum(abs, s∂y .- sfy) / eps() < 50

        # medium
        mdims = (21, 33)
        mvals = getpoints(f, Float64, mdims, FiniteDiffPlan)

        m∂x = zeros(Float64, length(mvals))
        m∂y = zeros(Float64, length(mvals))

        differentiate!(m∂x, m∂y, mvals, mdims)

        mfx = getpoints(fx, Float64, mdims, FiniteDiffPlan)
        @test maximum(abs, s∂x .- sfx) / eps() < 500

        mfy = getpoints(fy, Float64, mdims, FiniteDiffPlan)
        @test maximum(abs, s∂y .- sfy) / eps() < 500

        # large
        ldims = (331, 123)
        lvals = getpoints(f, Float64, ldims, FiniteDiffPlan)

        l∂x = zeros(Float64, length(lvals))
        l∂y = zeros(Float64, length(lvals))

        differentiate!(l∂x, l∂y, lvals, ldims)

        lfx = getpoints(fx, Float64, ldims, FiniteDiffPlan)
        @test maximum(abs, s∂x .- sfx) / eps() < 5000

        lfy = getpoints(fy, Float64, ldims, FiniteDiffPlan)
        @test maximum(abs, s∂x .- sfx) / eps() < 5000
    end

    let f(x, y) = sinpi(x) * exp(y)
        dims = (1321, 1235)
        vals = getpoints(f, Float64, dims, FiniteDiffPlan)

        ∂x = zeros(Float64, length(vals))
        ∂y = zeros(Float64, length(vals))

        differentiate!(∂x, ∂y, vals, dims)

        test∂x = getpoints(Float64, dims, FiniteDiffPlan) do x, y
            π * cospi(x) * exp(y)
        end

        test∂y = getpoints(Float64, dims, FiniteDiffPlan) do x, y
            sinpi(x) * exp(y)
        end

        @test maximum(abs, test∂x .- ∂x) < 1e-4
        @test maximum(abs, test∂y .- ∂y) < 1e-4
    end

    let f(x, y) = (x * exp(-x^2 - y^2), x + y, y * sinpi(x), x * y)
        dims = (75, 125)
        vals = getpoints(f, Float64, dims, FiniteDiffPlan)
        N = dims[1] * dims[2]

        ∂x = zeros(Float64, N, 4)
        ∂y = zeros(Float64, N, 4)

        differentiate!(eachcol(∂x), eachcol(∂y), vals, dims)

        for i in 1:4
            ∂xi = zeros(Float64, N)
            ∂yi = zeros(Float64, N)
            differentiate!(∂xi, ∂yi, vals[i], dims)

            @test ∂x[:, i] == ∂xi
            @test ∂y[:, i] == ∂yi
        end
    end
end

@safetestset "Simpson weights" begin
    using ..FiniteDifferences: getpoints, FiniteDiffPlan, getsimpweights

    @test getsimpweights(Float64, 5, 3) ≈ [
        1, 4, 1,  4, 16, 4,  2, 8, 2,  4, 16, 4,  1, 4, 1
    ] ./ (9 * 4 * 2)

    @test getsimpweights(Float64, 5, 7) ≈ [
        1  4  2  4  1
        4 16  8 16  4
        2  8  4  8  2
        4 16  8 16  4
        2  8  4  8  2
        4 16  8 16  4
        1  4  2  4  1
    ] ./ (9 * 4 * 6) |> vec

    let f(x, y) = 1 + 2x + 3y + 4x^2 + 5x*y + 6y^2
        nx, ny = (3, 3)
        vals = getpoints(f, Float64, (nx, ny), FiniteDiffPlan)
        weights = getsimpweights(Float64, nx, ny)
        @test sum(weights .* vals) ≈ 97 / 12 atol=10eps()
    end
end

@safetestset "compute!" begin
    using ..FiniteDifferences: getpoints, FiniteDiffPlan, getsimpweights
    using PoincareInvariants
    using LinearAlgebra: dot

    f(x, y) = [
         1 +  2*x +  3*y +  4*x^2 +  5*x*y +  6*y^2,
         7 +  8*x +  9*y + 10*x^2 + 11*x*y + 12*y^2,
        13 + 14*x + 15*y + 16*x^2 + 17*x*y + 18*y^2,
        19 + 20*x + 21*y + 22*x^2 + 23*x*y + 24*y^2]

    fx(x, y) = [
         2 +  8*x +  5*y,
         8 + 20*x + 11*y,
        14 + 32*x + 17*y,
        20 + 44*x + 23*y]

    fy(x, y) = [
         3 +  5*x + 12*y,
         9 + 11*x + 24*y,
        15 + 17*x + 36*y,
        21 + 23*x + 48*y]

    Ω(z, ::Any, ::Any) = [0 0 z[1] z[2]; 0 0 z[3] z[4]; -z[1] -z[3] 0 0; -z[2] -z[4] 0 0]
    dims = (11, 17)
    integrand = getpoints(Float64, dims, FiniteDiffPlan) do x, y
        dot(fy(x, y), Ω(f(x, y), 0, nothing), fx(x, y))
    end

    testI = dot(getsimpweights(Float64, dims...), integrand)

    pinv = SecondPoincareInvariant{Float64, 4}(Ω, (11, 17), FiniteDiffPlan)

    @test pinv.plan isa FiniteDiffPlan

    points = getpoints(f, Float64, dims, FiniteDiffPlan)
    I = compute!(pinv, points, 0, nothing)

    for i in 1:4
        @test pinv.plan.∂x[i] ≈ getpoints(fx, Float64, dims, FiniteDiffPlan)[i]
        @test pinv.plan.∂y[i] ≈ getpoints(fy, Float64, dims, FiniteDiffPlan)[i]
    end

    @test I ≈ testI atol=5e-11
end
