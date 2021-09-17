@safetestset "Unit Tests: Coefficient and Point Counts" begin
    using PoincareInvariants
    using PoincareInvariants.SecondPoincareInvariants: get_padua_num, get_n, check_padua_num,
        get_uu_coeff_num, get_uu_point_num

    for n in 1:11:10_000
        padua_num = get_padua_num(n)
        @test padua_num == (n + 1) * (n + 2) ÷ 2
        @test isinteger(get_n(padua_num))
        @test get_n(padua_num) == n
        @test (check_padua_num(padua_num); true)
        @test_throws ArgumentError check_padua_num(padua_num + 1)
        @test_throws ArgumentError check_padua_num(padua_num - 1)

        @test next_padua_num(padua_num) == padua_num
        @test next_padua_num(padua_num - 1) == padua_num
        @test next_padua_num(padua_num + 1) == get_padua_num(n + 1)
        
        @test get_uu_coeff_num(n) == n * (n + 1) ÷ 2
        @test get_uu_point_num(n) == n^2
        @test get_uu_coeff_num(n) == (get_uu_point_num(n) - n) ÷ 2 + n
    end
end

@safetestset "Unit Tests: 0..1 × 0..1 Square with degree 2" begin
    using PoincareInvariants
    using StaticArrays: SVector

    # We will use the following parameterisation as a test:
    # let f: (u, v) -> ((u + 1) / 2, (v + 1) / 2)
    # This maps the Padua points on -1..1 × -1..1 to points on 0..1 × 0..1
    # In other words, it maps the true Padua points onto our padua points on 0..1 × 0..1,
    # since get_padua_points(n) returns n points on 0..1 × 0..1
    #
    # For such a simple map we shall use 6 Padua points. These map onto the coefficient matrix
    # [f00 f01 f02;
    #  f10 f11  0 ;
    #  f20  0   0 ]
    # Where fij refers to Ti(u) * Tj(v) with Tn(x), the nth Chebyshev polynomial
    # The first three Chebyshev polynomials are
    # T0(x) = 1
    # T1(x) = x
    # T2(x) = 2x² - 1
    # Since our parameterisation is linear, going up to T1 would have been enough, too
    
    padua_num = 6
    paduapoints = get_padua_points(padua_num)

    # For testing we will add some extra dimensions, but the extra components will be set to 0
    @testset "SecondPoincareInvariant{$N, $T}($padua_num)" for T in [Float32, Float64], N in [2, 6]
        pinv = SecondPoincareInvariant{N, T}(padua_num)

        @test pinv isa AbstractPoincareInvariant

        @test pinv.degree == 2

        Ω = CanonicalSymplecticMatrix(N)

        idx = N ÷ 2 + 1

        phasepoints = zeros(T, padua_num, N)
        phasepoints[:, 1] .= paduapoints[:, 1]
        phasepoints[:, idx] .= paduapoints[:, 2]

        # The area of a 0..1 × 0..1 square should be 1. Is it?
        out = compute(pinv, phasepoints, Ω)
        @test out ≈ 1 atol=1e-10
        @test out isa T

        # Let's check our working
        # The matrix of coefficients gets turned into a vector as follows
        # [f00, f01, f10, f02, f11, f20]

        @test eltype(pinv.cc_coeffs) == T
        @test size(pinv.cc_coeffs) == (6, N)

        # For (u, v) -> x (first dimension), we have: x = 0.5 T00(u, v) + 0.5 T10(u, v)
        @test pinv.cc_coeffs[:, 1] ≈ [0.5, 0, 0.5, 0, 0, 0] atol=1e-10

        # For (u, v) -> y (second dimension), we have: y = 0.5 T00(u, v) + 0.5 T01(u, v)
        @test pinv.cc_coeffs[:, idx] ≈ [0.5, 0.5, 0, 0, 0, 0] atol=1e-10

        if N > 2
            # For all higher dimes, we have: z = 0
            @test all(isapprox.(pinv.cc_coeffs[:, 2:idx-1], 0, atol=1e-10))
            @test all(isapprox.(pinv.cc_coeffs[:, idx+1:end], 0, atol=1e-10))
        end

        # pinv.D1toUU and pinv.D1toUU differentiate and
        # convert to a basis using Chebyshev polynomials of the second kind
        # These are:
        # U0(x) = 1
        # U1(x) = 2x
        # U2(x) = 4x² - 1

        @test pinv.D1toUU == [0 0 1 0  0  0; # U00
                              0 0 0 0 0.5 0; # U01
                              0 0 0 0  0  2] # U10

        @test pinv.D2toUU == [0 1 0 0  0  0; # U00
                              0 0 0 2  0  0; # U01
                              0 0 0 0 0.5 0] # U10
        
        # pinv.CCtoUU converts from Chebyshev of first kind to second kind
        # Also truncates highest order coeffs
        @test pinv.CCtoUU == [1  0   0  -0.5 0 -0.5; # U00
                              0 0.5  0    0  0   0 ; # U01
                              0  0  0.5   0  0   0 ] # U10
        
        @test eltype(pinv.uu_coeffs) == T
        @test size(pinv.uu_coeffs) == (3, N)

        # x = 0.5 U00(u, v) + 0.25 U10(u, v)
        @test pinv.uu_coeffs[:, 1] ≈ [0.5, 0, 0.25] atol=1e-10

        # y = 0.5 U00(u, v) + 0.25 U01(u, v)
        @test pinv.uu_coeffs[:, idx] ≈ [0.5, 0.25, 0] atol=1e-10

        if N > 2
            # z = 0
            @test all(isapprox.(pinv.uu_coeffs[:, 2:idx-1], 0, atol=1e-10))
            @test all(isapprox.(pinv.uu_coeffs[:, idx+1:end], 0, atol=1e-10))
        end

        @test eltype(pinv.uu_d1_coeffs) == T
        @test size(pinv.uu_d1_coeffs) == (3, N)

        @test pinv.uu_d1_coeffs[:, 1] ≈ [0.5, 0, 0] atol=1e-10
        @test pinv.uu_d1_coeffs[:, idx] ≈ [ 0 , 0, 0] atol=1e-10

        @test eltype(pinv.uu_d2_coeffs) == T
        @test size(pinv.uu_d2_coeffs) == (3, N)

        @test pinv.uu_d2_coeffs[:, 1] ≈ [ 0 , 0, 0] atol=1e-10
        @test pinv.uu_d2_coeffs[:, idx] ≈ [0.5, 0, 0] atol=1e-10

        @test eltype(pinv.uu_points) == SVector{2, T}
        @test size(pinv.uu_points) == (4,)

        @test eltype(pinv.uu_vals) == T
        @test size(pinv.uu_vals) == (4, N)

        @test pinv.uu_vals[:, 1] ≈ [(v[1] + 1) / 2 for v in pinv.uu_points] atol=1e-10
        @test pinv.uu_vals[:, idx] ≈ [(v[2] + 1) / 2 for v in pinv.uu_points] atol=1e-10

        if N > 2
            @test all(isapprox.(pinv.uu_vals[:, 2:idx-1], 0, atol=1e-10))
            @test all(isapprox.(pinv.uu_vals[:, idx+1:end], 0, atol=1e-10))
        end

        @test eltype(pinv.uu_d1_vals) == T
        @test size(pinv.uu_d1_vals) == (4, N)

        @test pinv.uu_d1_vals[:, 1] ≈ [0.5, 0.5, 0.5, 0.5] atol=1e-10
        @test pinv.uu_d1_vals[:, idx] ≈ [0, 0, 0, 0] atol=1e-10

        if N > 2
            @test all(isapprox.(pinv.uu_d1_vals[:, 2:idx-1], 0, atol=1e-10))
            @test all(isapprox.(pinv.uu_d1_vals[:, idx+1:end], 0, atol=1e-10))
        end

        @test eltype(pinv.uu_d2_vals) == T
        @test size(pinv.uu_d2_vals) == (4, N)

        @test pinv.uu_d2_vals[:, 1] ≈ [0, 0, 0, 0] atol=1e-10
        @test pinv.uu_d2_vals[:, idx] ≈ [0.5, 0.5, 0.5, 0.5] atol=1e-10

        if N > 2
            @test all(isapprox.(pinv.uu_d2_vals[:, 2:idx-1], 0, atol=1e-10))
            @test all(isapprox.(pinv.uu_d2_vals[:, idx+1:end], 0, atol=1e-10))
        end

        @test eltype(pinv.uu_I_vals) == T
        @test pinv.uu_I_vals ≈ [0.25, 0.25, 0.25, 0.25] atol=1e-10

        @test eltype(pinv.uu_I_coeffs) == T
        @test pinv.uu_I_coeffs ≈ [0.25, 0, 0] atol=1e-10

        @test eltype(pinv.UUIntegral) == T
        @test pinv.UUIntegral == [4, 0, 0]'
    end
end

@safetestset "Function Tests: Free Particle" begin
    using PoincareInvariants

    function free_particle!(state, δt)
        mid = length(state) ÷ 2
        state[1:mid] += state[mid+1:end] .* δt
    end

    padua_num = next_padua_num(50 * 50)
    paduapoints = get_padua_points(padua_num)

    @testset "Free Particle in $(N)D on 0..1 × 0..1" for N in [2, 6]
        pinv = SecondPoincareInvariant{N, Float64}(padua_num)
        Ω = CanonicalSymplecticMatrix(N)

        phasepoints = ones(Float64, padua_num, N)
        phasepoints[:, 1] .= paduapoints[:, 1]
        phasepoints[:, N ÷ 2 + 1] .= paduapoints[:, 2]

        @test compute(pinv, phasepoints, Ω) ≈ 1 atol=1e-10

        free_particle!.(eachrow(phasepoints), 10)
        @test compute(pinv, phasepoints, Ω) ≈ 1 atol=1e-10

        free_particle!.(eachrow(phasepoints), 100)
        @test compute(pinv, phasepoints, Ω) ≈ 1 atol=1e-10

        free_particle!.(eachrow(phasepoints), 1000)
        @test compute(pinv, phasepoints, Ω) ≈ 1 atol=1e-10
    end

    @testset "Free Particle in $(N)D on Quarter Circle" for N in [2, 6]
        pinv = SecondPoincareInvariant{N, Float64}(padua_num)
        Ω = CanonicalSymplecticMatrix(N)

        phasepoints = ones(Float64, padua_num, N)
        phasepoints[:, 1] .= map(eachrow(paduapoints)) do row
            row[1] * cos(row[2] * π/2)
        end
        phasepoints[:, N ÷ 2 + 1] .= map(eachrow(paduapoints)) do row
            row[1] * sin(row[2] * π/2)
        end

        @test compute(pinv, phasepoints, Ω) ≈ π/4 atol=1e-10

        free_particle!.(eachrow(phasepoints), 10)
        @test compute(pinv, phasepoints, Ω) ≈ π/4 atol=1e-10

        free_particle!.(eachrow(phasepoints), 100)
        @test compute(pinv, phasepoints, Ω) ≈ π/4 atol=1e-10

        free_particle!.(eachrow(phasepoints), 1000)
        @test compute(pinv, phasepoints, Ω) ≈ π/4 atol=1e-10
    end
end

# @safetestset "Simple Surfaces" begin
#     using PoincareInvariants
#     using RandomMatrices
#     using LinearAlgebra: I

#     @testset "0..1 × 0..1 Square in $(N)D" for N in [2, 4, 10]
#         padua_num = next_padua_num(6)
#         pinv = SecondPoincareInvariant{N, Float64}(padua_num)
#         Ω = CanonicalSymplecticMatrix(N)

#         paduapoints = get_padua_points(padua_num)

#         phasepoints = zeros(Float64, padua_num, N)
#         phasepoints[:, 1:2] = paduapoints

#         @test compute(pinv, phasepoints, Ω) ≈ 1 rtol=1e-12

#         # 1 : 1 * 1
#         # 2 : 1 * y
#         # 3 : x * 1
#         @test pinv.cc_coeffs[:, 1] ≈ Float64[i == 1 || i == 3 ? 0.5 : 0 for i in 1:padua_num] atol=1e-10
#         @test pinv.cc_coeffs[:, 2] ≈ Float64[i == 1 || i == 2 ? 0.5 : 0 for i in 1:padua_num] atol=1e-10

#         N > 2 && (@test all(isapprox.(pinv.cc_coeffs[:, 3:end], 0, atol=1e-10)))

#         uu_coeff_num = size(pinv.uu_coeffs, 1)
#         # 1 : 1 * 1
#         # 2 : 1 * 2y
#         # 3 : 2x * 1
#         @test pinv.uu_coeffs[:, 1] ≈ map(1:uu_coeff_num) do i
#             i == 1 && return 0.5
#             i == 3 && return 0.25
#             return 0.0
#         end atol=1e-10
#         @test pinv.uu_coeffs[:, 2] ≈ map(1:uu_coeff_num) do i
#             i == 1 && return 0.5
#             i == 2 && return 0.25
#             return 0.0
#         end atol=1e-10

#         N > 2 && (@test all(isapprox.(pinv.uu_coeffs[:, 3:end], 0, atol=1e-10)))

#         @test pinv.D1toUU ≈ [0 0 0.5 0 0 0;
#                              0 0 0 0 0.5 0;
#                              0 0 0 0 0 2  ]

#         display(pinv.D1toUU)

#         @test pinv.uu_d1_coeffs[:, 1] ≈ Float64[i == 1 ? 2 * 0.25 : 0 for i in 1:uu_coeff_num] atol=1e-10
#         @test pinv.uu_d2_coeffs[:, 2] ≈ Float64[i == 1 ? 2 * 0.25 : 0 for i in 1:uu_coeff_num] atol=1e-10

#         @test pinv.uu_d1_coeffs[:, 2] ≈ zeros(uu_coeff_num) atol=1e-10
#         @test pinv.uu_d2_coeffs[:, 1] ≈ zeros(uu_coeff_num) atol=1e-10
#     end

#     # @testset "Plane in $(N)D" for N in [2, 4, 10]
#     #     padua_num = next_padua_num(10_000)
#     #     pinv = SecondPoincareInvariant{N, Float64}(padua_num)
#     #     Ω = CanonicalSymplecticMatrix(N)

#     #     # Pick random orthogonal matrix
#     #     # and select two vectors
#     #     Q = rand(Haar(1), N)
#     #     @assert Q * Q' ≈ I
#     #     P = Q[1:2, 1:N]

#     #     paduapoints = get_padua_points(padua_num)
        
#     #     phasepoints = paduapoints * P
#     #     @assert size(phasepoints) == (padua_num, N)

#     #     @test compute(pinv, phasepoints, Ω) ≈ 1 rtol=1e-12
#     # end
# end

# @safetestset "Free Particles" begin
#     using PoincareInvariants

#     @testset "$np Free Particles in 3D" for np in [2, 5], T in [Float32, Float64]
#         padua_num = next_padua_num(10_000)

#         N = np * 6

#         pinv = SecondPoincareInvariant{N, T}(padua_num)

#         function free_particle!(init, δt)
#             mid = length(init) ÷ 2
#             init[1:mid] += init[mid+1:end] .* δt
#         end

#         Ω = CanonicalSymplecticMatrix(N)

#         ppoints = get_padua_points(padua_num)
#         phasepoints = ones(T, padua_num, N)

#         for i in 1:padua_num
#             phasepoints[i, np * 3 .+ (1:2)] .+= ppoints[i, :]
#         end

#         @test compute(pinv, phasepoints, Ω) ≈ 1 rtol=1e-12

#         free_particle!.(eachrow(phasepoints), 10)
#         @test compute(pinv, phasepoints, Ω) ≈ 1 rtol=1e-12

#         free_particle!.(eachrow(phasepoints), 10)
#         @test compute(pinv, phasepoints, Ω) ≈ 1 rtol=1e-12
#     end
# end