using FastTransforms: PaduaTransformPlan

@testset "Coefficient and Point Counts" begin
    for n in 1:11:10_000
        padua_num = PoincareInvariants.get_padua_num(n)
        @test padua_num == (n + 1) * (n + 2) ÷ 2
        @test isinteger(PoincareInvariants.get_n(padua_num))
        @test PoincareInvariants.get_n(padua_num) == n
        @test (PoincareInvariants.check_padua_num(padua_num); true)
        @test_throws ArgumentError PoincareInvariants.check_padua_num(padua_num + 1)
        @test_throws ArgumentError PoincareInvariants.check_padua_num(padua_num - 1)

        @test next_padua_num(padua_num) == padua_num
        @test next_padua_num(padua_num - 1) == padua_num
        @test next_padua_num(padua_num + 1) == PoincareInvariants.get_padua_num(n + 1)
        
        @test PoincareInvariants.get_uu_coeff_num(n) == n * (n + 1) ÷ 2
        @test PoincareInvariants.get_uu_point_num(n) == n^2
        @test PoincareInvariants.get_uu_coeff_num(n) == (PoincareInvariants.get_uu_point_num(n) - n) ÷ 2 + n
    end
end

@testset "PoincareInvariant2 Struct" begin
    @test PoincareInvariant2 <: AbstractPoincareInvariant

    nums = next_padua_num.([10_000, 20_000])

    @testset "PoincareInvariant2{$N, $T}($padua_num)" for T in [Float32, Float64], N in [2, 36], padua_num in nums
        n = PoincareInvariants.get_n(padua_num)
        uu_coeff_num = PoincareInvariants.get_uu_coeff_num(n)
        uu_point_num = PoincareInvariants.get_uu_point_num(n)

        pinv = PoincareInvariant2{N, T}(padua_num)

        @test pinv isa PoincareInvariant2

        @test pinv.n == n

        @test pinv.cc_coeffs isa Matrix{T}
        @test size(pinv.cc_coeffs) == (padua_num, N)

        let v = [cos(v[1]) * sin(v[2]) for v in eachrow(get_padua_points(padua_num))]
            coeffs = pinv.padua_plan * copy(v)
            @test ipaduatransform(coeffs, Val{false}) ≈ v rtol=1e-6
        end

        @test size(pinv.D1toUU) == (uu_coeff_num, padua_num)
        @test size(pinv.D2toUU) == (uu_coeff_num, padua_num)
        @test size(pinv.CCtoUU) == (uu_coeff_num, padua_num)

        @test size(pinv.uu_coeffs) == (uu_coeff_num, N)
        @test size(pinv.uu_d1_coeffs) == (uu_coeff_num, N)
        @test size(pinv.uu_d2_coeffs) == (uu_coeff_num, N)

        @test eltype(pinv.uu_coeffs) == T
        @test eltype(pinv.uu_d1_coeffs) == T
        @test eltype(pinv.uu_d2_coeffs) == T

        @test size(pinv.uu_points) == (uu_point_num,)
        @test eltype(pinv.uu_points) == SVector{2, T}

        @test size(pinv.uu_vals) == (uu_point_num, N)
        @test size(pinv.uu_d1_vals) == (uu_point_num, N)
        @test size(pinv.uu_d2_vals) == (uu_point_num, N)

        @test eltype(pinv.uu_vals) == T
        @test eltype(pinv.uu_d1_vals) == T
        @test eltype(pinv.uu_d2_vals) == T

        let v = [cos(v[1]) * sin(v[2]) for v in pinv.uu_points]
            coeffs = pinv.uu_plan * copy(v)
            @test pinv.uu_iplan * copy(coeffs) ≈ v rtol=1e-6
        end

        @test size(pinv.uu_I_vals) == (uu_point_num,)
        @test eltype(pinv.uu_I_vals) == T

        @test size(pinv.uu_I_coeffs) == (uu_coeff_num,)
        @test eltype(pinv.uu_I_coeffs) == T

        @test size(pinv.UUIntegral) == (1, uu_coeff_num)
        @test eltype(pinv.UUIntegral) == T

        let v = [v[1] + 4 * v[2]^3 for v in pinv.uu_points]
            coeffs = pinv.uu_plan * copy(v)
            @test dot(pinv.UUIntegral, coeffs) ≈ 0.5 + 1 rtol=1e-6
        end
    end
end

@testset "Free Particles" begin
    @testset "$np $(nd)D Free Particle(s)" for np in 1:3, nd in 1:3
        padua_num = next_padua_num(10_000)
        ppoints = get_padua_points(padua_num)

        pinv = PoincareInvariant2{2, Float64}(padua_num)

        function free_particle!(init, δt)
            init[1] += init[2] .* δt
        end

        Ω = [0 -1; 1 0]

        @test compute(pinv, ppoints, Ω) ≈ 1

        free_particle!.(eachrow(ppoints), 10)
        @test compute(pinv, ppoints, Ω) ≈ 1

        free_particle!.(eachrow(ppoints), 10)
        @test compute(pinv, ppoints, Ω) ≈ 1
    end
end