using FastTransforms: paduapoints, plan_paduatransform!

@testset "Padua Transforms" begin
    @testset "get_padua_points" begin
        for n in 1:100
            N = (n + 1) * (n + 2) ÷ 2
            points = map(eachrow(paduapoints(n))) do row
                (SVector{2, Float64}(row) .+ 1) ./ 2
            end

            @test get_padua_points(N) == points
        end
    end

    @testset "paduatransform!" begin
        f(v) = v[1] * sin(v[2])

        padua_points = get_padua_points(105)
        v = f.(padua_points)

        plan = plan_paduatransform!(v, Val{false})

        coeffs = plan * copy(v)

        out = similar(v)
        PoincareInvariants.padua_transform!(out, v, plan)

        @test out == coeffs
        @test v == f.(padua_points)
    end
end  # big testset