@safetestset "CanoincalSymplecticMatrices" begin include("test_CanonicalSymplecticMatrices.jl") end

@safetestset "Chebyshev" begin include("test_Chebyshev.jl") end

@safetestset "Unit Tests for Implementation 1" begin
    using PoincareInvariants
    using PoincareInvariants.SecondPoincareInvariants: getpaduapoints, Setup1
    using StaticArrays: SVector

    # We will use the following parameterisation as a test:
    # let f: (u, v) -> ((u + 1) / 2, (v + 1) / 2)
    # This maps the Padua points on -1..1 × -1..1 to points on 0..1 × 0..1
    # In other words, it maps the true Padua points onto our Padua points on 0..1 × 0..1,
    # since getpaduapoints returns points on 0..1 × 0..1
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
    
    degree = 2
    paduanum = 6
    paduapoints = getpaduapoints(2)

    # For testing we will add some extra dimensions, but the extra components will be set to 0
    @testset "Setup1($D, $T, Ω, $paduanum)" for D in [2, 6], T in [Float32, Float64]
        Ω = CanonicalSymplecticMatrix(D)
        setup = Setup1(D, T, Ω, paduanum)
        pinv = SecondPoincareInvariant{D, T}(Ω, paduanum, setup)

        @test pinv isa AbstractPoincareInvariant

        # @test pinv.degree == 2

        # index of first momentum dimension
        idx = D ÷ 2 + 1

        phasepoints = zeros(T, paduanum, D)
        phasepoints[:, 1] .= paduapoints[:, 1]
        phasepoints[:, idx] .= paduapoints[:, 2]

        # The area of a 0..1 × 0..1 square should be 1. Is it?
        out = compute(pinv, phasepoints)
        @test out ≈ 1 atol=1e-10
        @test out isa T

        # Let's check our working
        # The matrix of coefficients gets turned into a vector as follows
        # [f00, f01, f10, f02, f11, f20]

        @test eltype(setup.C1coeffs) == T
        @test size(setup.C1coeffs) == (paduanum, D)

        # For (u, v) -> x (first dimension), we have: x = 0.5 T00(u, v) + 0.5 T10(u, v)
        @test setup.C1coeffs[:, 1] ≈ [0.5, 0, 0.5, 0, 0, 0] atol=1e-10

        # For (u, v) -> y (second dimension), we have: y = 0.5 T00(u, v) + 0.5 T01(u, v)
        @test setup.C1coeffs[:, idx] ≈ [0.5, 0.5, 0, 0, 0, 0] atol=1e-10

        if D > 2
            # For all higher dimes = 0
            @test all(isapprox.(setup.C1coeffs[:, 2:idx-1], 0, atol=1e-10))
            @test all(isapprox.(setup.C1coeffs[:, idx+1:end], 0, atol=1e-10))
        end

        # setup.DutoC2 and setup.DvtoC2 differentiate and
        # convert to a basis using Chebyshev polynomials of the second kind
        # These are:
        # U0(x) = 1
        # U1(x) = 2x
        # U2(x) = 4x² - 1

        @test setup.DutoC2 == [0 0 1 0  0  0; # U00
                               0 0 0 0 0.5 0; # U01
                               0 0 0 0  0  2] # U10

        @test setup.DvtoC2 == [0 1 0 0  0  0; # U00
                               0 0 0 2  0  0; # U01
                               0 0 0 0 0.5 0] # U10
        
        # setup.C1toC2 converts from Chebyshev of first kind to second kind
        # Also truncates highest order coeffs
        @test setup.C1toC2 == [1  0   0  -0.5 0 -0.5; # U00
                               0 0.5  0    0  0   0 ; # U01
                               0  0  0.5   0  0   0 ] # U10
        
        @test eltype(setup.C2coeffs) == T
        @test size(setup.C2coeffs) == (3, D)

        # x = 0.5 U00(u, v) + 0.25 U10(u, v)
        @test setup.C2coeffs[:, 1] ≈ [0.5, 0, 0.25] atol=1e-10

        # y = 0.5 U00(u, v) + 0.25 U01(u, v)
        @test setup.C2coeffs[:, idx] ≈ [0.5, 0.25, 0] atol=1e-10

        if D > 2
            # z = 0
            @test all(isapprox.(setup.C2coeffs[:, 2:idx-1], 0, atol=1e-10))
            @test all(isapprox.(setup.C2coeffs[:, idx+1:end], 0, atol=1e-10))
        end

        @test eltype(setup.C2Ducoeffs) == T
        @test size(setup.C2Ducoeffs) == (3, D)

        @test setup.C2Ducoeffs[:, 1] ≈ [0.5, 0, 0] atol=1e-10
        @test setup.C2Ducoeffs[:, idx] ≈ [ 0 , 0, 0] atol=1e-10

        @test eltype(setup.C2Dvcoeffs) == T
        @test size(setup.C2Dvcoeffs) == (3, D)

        @test setup.C2Dvcoeffs[:, 1] ≈ [ 0 , 0, 0] atol=1e-10
        @test setup.C2Dvcoeffs[:, idx] ≈ [0.5, 0, 0] atol=1e-10

        @test eltype(setup.C2vals) == T
        @test size(setup.C2vals) == (4, D)

        # @test setup.C2vals[:, 1] ≈ [(v[1] + 1) / 2 for v in pinv.uu_points] atol=1e-10
        # @test setup.C2vals[:, idx] ≈ [(v[2] + 1) / 2 for v in pinv.uu_points] atol=1e-10

        if D > 2
            @test all(isapprox.(setup.C2vals[:, 2:idx-1], 0, atol=1e-10))
            @test all(isapprox.(setup.C2vals[:, idx+1:end], 0, atol=1e-10))
        end

        @test eltype(setup.C2Duvals) == T
        @test size(setup.C2Duvals) == (4, D)

        @test setup.C2Duvals[:, 1] ≈ [0.5, 0.5, 0.5, 0.5] atol=1e-10
        @test setup.C2Duvals[:, idx] ≈ [0, 0, 0, 0] atol=1e-10

        if D > 2
            @test all(isapprox.(setup.C2Duvals[:, 2:idx-1], 0, atol=1e-10))
            @test all(isapprox.(setup.C2Duvals[:, idx+1:end], 0, atol=1e-10))
        end

        @test eltype(setup.C2Dvvals) == T
        @test size(setup.C2Dvvals) == (4, D)

        @test setup.C2Dvvals[:, 1] ≈ [0, 0, 0, 0] atol=1e-10
        @test setup.C2Dvvals[:, idx] ≈ [0.5, 0.5, 0.5, 0.5] atol=1e-10

        if D > 2
            @test all(isapprox.(setup.C2Dvvals[:, 2:idx-1], 0, atol=1e-10))
            @test all(isapprox.(setup.C2Dvvals[:, idx+1:end], 0, atol=1e-10))
        end

        @test eltype(setup.C2Ivals) == T
        @test setup.C2Ivals ≈ [0.25, 0.25, 0.25, 0.25] atol=1e-10

        @test eltype(setup.C2Icoeffs) == T
        @test setup.C2Icoeffs ≈ [0.25, 0, 0] atol=1e-10

        @test eltype(setup.C2Integral) == T
        @test setup.C2Integral == [4, 0, 0]'
    end
end

@safetestset "Free Particle" begin
    using PoincareInvariants

    function free_particle!(state, δt)
        mid = length(state) ÷ 2
        state[1:mid] += state[mid+1:end] .* δt
    end

    @testset "Free Particle in $(D)D on 0..1 × 0..1" for D in [2, 6]
        Ω = CanonicalSymplecticMatrix(D)
        pinv = SecondPoincareInvariant{D, Float64}(Ω, 500)

        parampoints = getpoints(pinv)

        phasepoints = ones(Float64, pinv.N, D)
        phasepoints[:, 1] .= parampoints[:, 1]
        phasepoints[:, D ÷ 2 + 1] .= parampoints[:, 2]

        @test compute(pinv, phasepoints) ≈ 1 atol=1e-10

        free_particle!.(eachrow(phasepoints), 10)
        @test compute(pinv, phasepoints) ≈ 1 atol=1e-10

        free_particle!.(eachrow(phasepoints), 100)
        @test compute(pinv, phasepoints) ≈ 1 atol=1e-10

        free_particle!.(eachrow(phasepoints), 1000)
        @test compute(pinv, phasepoints) ≈ 1 atol=1e-10
    end

    @testset "Free Particle in $(D)D on Quarter Circle" for D in [4, 12]
        Ω = CanonicalSymplecticMatrix(D)
        pinv = SecondPoincareInvariant{D, Float64}(Ω, 500)

        parampoints = getpoints(pinv)

        phasepoints = ones(Float64, pinv.N, D)
        phasepoints[:, 2] .= map(eachrow(parampoints)) do row
            row[1] * cos(row[2] * π/2)
        end
        phasepoints[:, D ÷ 2 + 2] .= map(eachrow(parampoints)) do row
            row[1] * sin(row[2] * π/2)
        end

        @test compute(pinv, phasepoints) ≈ π/4 atol=1e-10

        free_particle!.(eachrow(phasepoints), 10)
        @test compute(pinv, phasepoints) ≈ π/4 atol=1e-10

        free_particle!.(eachrow(phasepoints), 100)
        @test compute(pinv, phasepoints) ≈ π/4 atol=1e-10

        free_particle!.(eachrow(phasepoints), 1000)
        @test compute(pinv, phasepoints) ≈ π/4 atol=1e-10
    end
end
