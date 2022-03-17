using PoincareInvariants.SecondPoincareInvariants.FiniteDifferences

@safetestset "getpoints and getpointnum" begin
    using ..FiniteDifferences: FiniteDiffPlan, getpoints, getpointnum

    @test getpointnum(50, FiniteDiffPlan) == 64
    @test getpointnum(87, FiniteDiffPlan) == 100
    @test getpointnum(400, FiniteDiffPlan) == 400

    for N in [3123, 63287, 762384]
        @test getpointnum(N, FiniteDiffPlan) == ceil(Int, sqrt(N))^2
    end

    @test getpointnum((3, 5), FiniteDiffPlan) == 15
    @test getpointnum((423789, 326), FiniteDiffPlan) == 423789 * 326
    @test getpointnum((123, 908342), FiniteDiffPlan) == 123 * 908342

    @inferred NTuple{3, Vector{Float64}} getpoints((x, y) -> (x, y, x+y), Float64, (7, 4), FiniteDiffPlan)
    @inferred Vector{Float32} getpoints((x, y) -> sin(x), Float32, 64, FiniteDiffPlan)

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

    @test getpoints((x, y) -> y, Float64, 120, FiniteDiffPlan) ≈ Float64[
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

    testpnts = getpoints((x, y) -> (x, y), Float64, 6342, FiniteDiffPlan)
    pnts4 = getpoints(Float64, 6342, FiniteDiffPlan) do x, y
        x, y, x * y, x + y
    end
    @test pnts4 isa NTuple{4, Vector{Float64}}
    @test all(pnts4) do v
        length(v) == getpointnum(6342, FiniteDiffPlan)
    end

    @test pnts4[1] ≈ testpnts[1]
    @test pnts4[2] ≈ testpnts[2]
    @test pnts4[3] ≈ testpnts[1] .* testpnts[2]
    @test pnts4[4] ≈ testpnts[1] .+ testpnts[2]
end
