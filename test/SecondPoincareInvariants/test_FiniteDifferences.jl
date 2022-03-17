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
end
