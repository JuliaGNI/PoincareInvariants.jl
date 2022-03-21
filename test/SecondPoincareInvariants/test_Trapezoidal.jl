using PoincareInvariants.SecondPoincareInvariants.Trapezoidal

@safetestset "getpoints, getpointspec and getpointnum" begin
    using ..Trapezoidal: TrapezoidalPlan, getpoints, getpointnum, getpointspec

    @test getpointspec(50, TrapezoidalPlan) == (8, 8)
    @test getpointspec(87, TrapezoidalPlan) == (10, 10)
    @test getpointspec((53, 42), TrapezoidalPlan) == (53, 42)

    @test getpointnum(50, TrapezoidalPlan) == 64
    @test getpointnum(87, TrapezoidalPlan) == 100
    @test getpointnum(400, TrapezoidalPlan) == 400

    for N in [3123, 63287, 762384]
        @test getpointnum(N, TrapezoidalPlan) == ceil(Int, sqrt(N))^2
    end

    @test getpointnum((3, 5), TrapezoidalPlan) == 15
    @test getpointnum((423789, 326), TrapezoidalPlan) == 423789 * 326
    @test getpointnum((123, 908342), TrapezoidalPlan) == 123 * 908342

    @inferred NTuple{3, Vector{Float64}} getpoints((x, y) -> (x, y, x+y), Float64, (7, 4), TrapezoidalPlan)
    @inferred Vector{Float32} getpoints((x, y) -> sin(x), Float32, (8, 8), TrapezoidalPlan)

    @test getpoints((x, y) -> (x, y), Float64, (5, 3), TrapezoidalPlan)[1] ≈ [
        0   , 0   , 0   ,
        0.25, 0.25, 0.25,
        0.5 , 0.5 , 0.5 ,
        0.75, 0.75, 0.75,
        1.0 , 1.0 , 1.0 ]
    @test getpoints((x, y) -> (x, y), Float64, (5, 3), TrapezoidalPlan)[2] ≈ [
        0  , 0.5, 1.0,
        0  , 0.5, 1.0,
        0  , 0.5, 1.0,
        0  , 0.5, 1.0,
        0  , 0.5, 1.0]

    reshapepnts = getpoints((x, y) -> (3x, 4y), Float64, (4, 5), TrapezoidalPlan)
    @test reshape(reshapepnts[1], 5, 4) ≈ [
        0 1 2 3;
        0 1 2 3;
        0 1 2 3;
        0 1 2 3;
        0 1 2 3]
    @test reshape(reshapepnts[2], 5, 4) ≈ [
        0 0 0 0;
        1 1 1 1;
        2 2 2 2;
        3 3 3 3;
        4 4 4 4]

    @test getpoints((x, y) -> y, Float64, (11, 11), TrapezoidalPlan) ≈ Float64[
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

    ps = getpointspec(6342, TrapezoidalPlan)
    testpnts = getpoints((x, y) -> (x, y), Float64, ps, TrapezoidalPlan)
    pnts4 = getpoints(Float64, ps, TrapezoidalPlan) do x, y
        x, y, x * y, x + y
    end
    @test pnts4 isa NTuple{4, Vector{Float64}}
    @test all(pnts4) do v
        length(v) == getpointnum(ps, TrapezoidalPlan)
    end

    @test pnts4[1] ≈ testpnts[1]
    @test pnts4[2] ≈ testpnts[2]
    @test pnts4[3] ≈ testpnts[1] .* testpnts[2]
    @test pnts4[4] ≈ testpnts[1] .+ testpnts[2]
end

@safetestset "Differentiation" begin
    using ..Trapezoidal: getpoints, TrapezoidalPlan, differentiate!

    let f(x, y) = sinpi(x) * exp(-y^2) + x * y
        dims = (4, 3)
        vals = getpoints(f, Float64, dims, TrapezoidalPlan)

        ∂x = zeros(Float64, length(vals))
        ∂y = zeros(Float64, length(vals))

        differentiate!(∂x, ∂y, vals, dims)

        Δx = 1 / 3
        @test ∂x ≈ [
            f(1/3,  0)-f(0,  0)  (f(2/3,  0)-f(0,  0))/2  (f(1,  0)-f(1/3,  0))/2  f(1,  0)-f(2/3,  0);
            f(1/3,1/2)-f(0,1/2)  (f(2/3,1/2)-f(0,1/2))/2  (f(1,1/2)-f(1/3,1/2))/2  f(1,1/2)-f(2/3,1/2);
            f(1/3,  1)-f(0,  1)  (f(2/3,  1)-f(0,  1))/2  (f(1,  1)-f(1/3,  1))/2  f(1,  1)-f(2/3,  1);
        ] ./ Δx |> vec

        Δy = 1 / 2
        @test ∂y ≈ [
            f(0,1/2)-f(0,  0)     f(1/3,1/2)-f(1/3,  0)     f(2/3,1/2)-f(2/3,  0)     f(1,1/2)-f(1,  0);
            (f(0,  1)-f(0,  0))/2 (f(1/3,  1)-f(1/3,  0))/2 (f(2/3,  1)-f(2/3,  0))/2 (f(1,  1)-f(1,  0))/2;
            f(0,  1)-f(0,1/2)     f(1/3,  1)-f(1/3,1/2)     f(2/3,  1)-f(2/3,1/2)     f(1,  1)-f(1,1/2);
        ] ./ Δy |> vec
    end

    let f(x, y) = sinpi(x) * y^2
        dims = (1000, 1000)
        vals = getpoints(f, Float64, dims, TrapezoidalPlan)

        ∂x = zeros(Float64, length(vals))
        ∂y = zeros(Float64, length(vals))

        differentiate!(∂x, ∂y, vals, dims)

        test∂x = getpoints(Float64, dims, TrapezoidalPlan) do x, y
            π * cospi(x) * y^2
        end

        test∂y = getpoints(Float64, dims, TrapezoidalPlan) do x, y
            sinpi(x) * 2 * y
        end

        @test maximum(abs, test∂y .- ∂y) < 2e-3
        @test maximum(abs, test∂y .- ∂y) < 2e-3
    end

    let f(x, y) = (x * exp(-x^2 - y^2), x + y, y * sinpi(x), x * y)
        dims = (75, 125)
        vals = getpoints(f, Float64, dims, TrapezoidalPlan)
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
