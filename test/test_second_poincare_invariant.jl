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

@testset "SecondPoincareInvariant Struct" begin
    @test SecondPoincareInvariant <: AbstractPoincareInvariant

    nums = next_padua_num.([10_000, 20_000])

    @testset "SecondPoincareInvariant{$N, $T}($padua_num)" for T in [Float32, Float64], N in [2, 36], padua_num in nums
        n = PoincareInvariants.get_n(padua_num)
        uu_coeff_num = PoincareInvariants.get_uu_coeff_num(n)
        uu_point_num = PoincareInvariants.get_uu_point_num(n)

        pinv = SecondPoincareInvariant{N, T}(padua_num)

        @test pinv isa SecondPoincareInvariant

        @test pinv.n == n

        @test pinv.cc_coeffs isa Matrix{T}
        @test size(pinv.cc_coeffs) == (padua_num, N)

        let v = [cos(v[1]) * sin(v[2]) for v in eachrow(get_padua_points(padua_num))]
            coeffs = pinv.padua_plan * copy(v)
            @test ipaduatransform(coeffs, Val{false}) ≈ v rtol=√eps(T)
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
            @test pinv.uu_iplan * copy(coeffs) ≈ v rtol=√eps(T)
        end

        @test size(pinv.uu_I_vals) == (uu_point_num,)
        @test eltype(pinv.uu_I_vals) == T

        @test size(pinv.uu_I_coeffs) == (uu_coeff_num,)
        @test eltype(pinv.uu_I_coeffs) == T

        @test size(pinv.UUIntegral) == (1, uu_coeff_num)
        @test eltype(pinv.UUIntegral) == T

        let v = [v[1] + 4 * v[2]^3 for v in pinv.uu_points]
            coeffs = pinv.uu_plan * copy(v)
            @test dot(pinv.UUIntegral, coeffs) ≈ 0.5 + 1 rtol=√eps(T)
        end
    end
end

@testset "Free Particles" begin
    @testset "$np Free Particles in 3D" for np in [2, 5], T in [Float32, Float64]
        padua_num = next_padua_num(10_000)

        N = np * 6

        pinv = SecondPoincareInvariant{N, T}(padua_num)

        function free_particle!(init, δt)
            mid = length(init) ÷ 2
            init[1:mid] += init[mid+1:end] .* δt
        end

        Ω = BlockArray(zeros(Int, N, N), [N ÷ 2, N ÷ 2], [N ÷ 2, N ÷ 2])
        Ω[Block(1, 2)] - I
        Ω[Block(2, 1)] + I

        ppoints = get_padua_points(padua_num)
        phasepoints = ones(T, padua_num, N)

        for i in 1:padua_num
            phasepoints[i, end-1:end] .+= ppoints[i, 1:2]
        end

        @test compute(pinv, phasepoints, Ω) ≈ 1

        free_particle!.(eachrow(phasepoints), 10)
        @test compute(pinv, phasepoints, Ω) ≈ 1 rtol=1e-12

        free_particle!.(eachrow(phasepoints), 10)
        @test compute(pinv, phasepoints, Ω) ≈ 1
    end
end