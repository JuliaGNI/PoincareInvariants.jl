using Test, PoincareInvariants

using Random: MersenneTwister

using StaticArrays

@testset "PoincareInvariants.jl" begin
    @testset "get_phase_points!" begin
        point_num = 100
        T = Float64
        N = 4

        param_func(v) = SVector{4}(
            v[1] * v[2],
            v[1] * 2,
            v[1] + v[2],
            v[1] + 2
        )

        out = ntuple(_ -> zeros(T, point_num), 4)

        param_points = let rng = MersenneTwister(1234)
            [SVector(rand(rng, T), rand(rng, T)) for _ in 1:point_num]
        end

        PoincareInvariants.get_phase_points!(param_func, out, param_points, Val(4))

        @test out[begin] == first.(param_func.(param_points))

        @test out[end] == last.(param_func.(param_points))
    end

    @testset "PoincareInvariant2" begin
        @test PoincareInvariant2 <: AbstractPoincareInvariant

        # pick odd type to make sure types are correct
        param_func(v) = SVector{2, Float32}(v[1] * v[2], v[1] + v[2])
        min_point_num = 321

        pinv = PoincareInvariant2{2}(param_func, min_point_num)

        @test pinv isa PoincareInvariant2

        @test pinv.phase_points isa NTuple{2, Vector{Float32}}

        @test all(length.(pinv.phase_points) .== pinv.point_num .>= min_point_num)
    end
end
