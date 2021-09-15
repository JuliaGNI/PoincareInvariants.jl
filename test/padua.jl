@testset "Padua Transforms" begin
    @testset "get_padua_points" begin
        @test eltype(get_padua_points(21)) == Float64

        @testset "get_padua_points($T, $n)" for T in [Float32, Float64], n in 1:11:100
            N = (n + 1) * (n + 2) ÷ 2
            ft_points = (paduapoints(T, n) .+ 1) ./ 2
            pi_points = get_padua_points(T, N)
            
            @test pi_points ≈ ft_points
            @test eltype(pi_points) == T
        end
    end

    @testset "paduatransform!" begin
        f(v) = v[1] * sin(v[2])

        padua_points = get_padua_points(105)
        v = f.(eachrow(padua_points))

        plan = plan_paduatransform!(v, Val{false})

        coeffs = plan * copy(v)

        out = similar(v)
        PoincareInvariants.padua_transform!(out, v, plan)

        @test out == coeffs
        @test v == f.(eachrow(padua_points))
    end
end  # big testset