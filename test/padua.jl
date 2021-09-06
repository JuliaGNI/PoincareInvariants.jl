using FastTransforms: paduapoints, plan_paduatransform!

@testset "Padua Transforms" begin
    @testset "get_padua_points" begin
        for n in 1:5:100
            N = (n + 1) * (n + 2) ÷ 2
            points = SVector{2, Float64}.(eachrow(paduapoints(n)))

            @test get_padua_points(N) == points
            @test get_padua_points(N - 1) == points
        end
    end

    @testset "paduatransform!" begin
        f(v) = v[1] * sin(v[2])

        padua_points = get_padua_points(100)
        v = f.(padua_points)

        plan = plan_paduatransform!(v, Val{false})

        coeffs = plan * copy(v)

        out = similar(v)
        PoincareInvariants.padua_transform!(out, v, plan)

        @test out == coeffs
        @test v == f.(padua_points)
    end
end  # big testset