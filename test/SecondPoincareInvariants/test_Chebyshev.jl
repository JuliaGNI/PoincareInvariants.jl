@safetestset "PaduaTransforms" begin include("test_PaduaTransforms.jl") end

using PoincareInvariants.SecondPoincareInvariants.Chebyshev

@safetestset "Differentiation" begin
    using ..Chebyshev: DiffPlan, differentiate!

    @test DiffPlan{Float64}(5).D == [0  1  0  3  0   5;
                                     0  0  4  0  8   0;
                                     0  0  0  6  0  10;
                                     0  0  0  0  8   0;
                                     0  0  0  0  0  10;
                                     0  0  0  0  0   0]

    @test DiffPlan{Float64}(4).D == [0  1  0  3  0;
                                     0  0  4  0  8;
                                     0  0  0  6  0;
                                     0  0  0  0  8;
                                     0  0  0  0  0]

    coeffs = [1  2  3  0;
              4  5  0  0;
              6  0  0  0;
              0  0  0  0]

    P = DiffPlan{Float64}(3)

    ∂x = zeros(4, 4)
    ∂y = zeros(4, 4)

    differentiate!(∂x, ∂y, P, coeffs)

    @test ∂x == [2  12  0  0;
                 5   0  0  0;
                 0   0  0  0;
                 0   0  0  0]

    @test ∂y == [ 4  5  0  0;
                 24  0  0  0;
                  0  0  0  0;
                  0  0  0  0]

    ndcoeffs = ntuple(_ -> rand(4, 4), 6)
    nd∂x = ntuple(_ -> zeros(4, 4), 6)
    nd∂y = ntuple(_ -> zeros(4, 4), 6)

    differentiate!(nd∂x, nd∂y, P, ndcoeffs)

    for i in 1:6
        differentiate!(∂x, ∂y, P, ndcoeffs[i])
        @test nd∂x[i] == ∂x
        @test nd∂y[i] == ∂y
    end
end

@safetestset "Integration" begin
    using ..Chebyshev:
        getintweights, integrate

    # integrating odd polynomials over symmetric boundary conditions gives 0
    @test all(getintweights(100)[2:2:end] .== 0)
    @test all(getintweights(999)[2:2:end] .== 0)

    coeffs = [1 4 7 0;
              2 5 8 0;
              3 6 9 0;
		      0 0 0 0]

    @test integrate(coeffs, getintweights(3)) ≈ 1 * 4 + 3 * -4/3 + 7 * -4/3 + 9 * 4/9 atol=5eps()
end

@safetestset "getintegrand! with CallIntPlan" begin
    using ..Chebyshev.PaduaTransforms
    using ..Chebyshev:
        DiffPlan, differentiate!, CallIntPlan, getintegrand!
    using PoincareInvariants.CanonicalSymplecticStructures

    @testset "6 Points in 2 Dimensions" begin
        degree = 2
        D = 2
        phasecoeffs = (
            Float64[1 2 3;
                    4 5 0;
                    6 0 0],
            Float64[1 4 6;
                    2 5 0;
                    3 0 0]
        )

        ∂xcoeffs = ntuple(_ -> zeros(degree+1, degree+1), D)
        ∂ycoeffs = ntuple(_ -> zeros(degree+1, degree+1), D)
        differentiate!(∂xcoeffs, ∂ycoeffs, DiffPlan{Float64}(degree), phasecoeffs)

        phasepoints = [Vector{Float64}(undef, getpaduanum(degree)) for _ in 1:D]
        iplan = InvPaduaTransformPlan{Float64}(degree)
        invpaduatransform!(phasepoints, iplan, phasecoeffs)

        Ω(v, t, p) = [0 -1; 1  0]

        plan = CallIntPlan{Float64}(D, degree)

        intcoeffs = zeros(degree+1, degree+1)

        getintegrand!(intcoeffs, plan, Ω, phasepoints, 0, nothing, ∂xcoeffs, ∂ycoeffs)

        @test intcoeffs ≈ [-72   -82  -30;
                           -82  -432    0;
                           -30     0    0] atol=5eps()
    end

    @testset "50 Points in 6 Dimensions" begin
        D = 6
        N = 50
        pointnum = nextpaduanum(N)
        degree = getdegree(pointnum)

        phasecoeffs = ntuple(_ -> zeros(degree+1, degree+1), D)
        phasecoeffs[1][1, 1] = 0.5  # 0.5
        phasecoeffs[1][1, 2] = 0.5  # 0.5 x
        phasecoeffs[D ÷ 2 + 1][1, 1] = 0.5  # 0.5
        phasecoeffs[D ÷ 2 + 1][2, 1] = 0.5  # 0.5 y

        ∂xcoeffs = ntuple(_ -> zeros(degree+1, degree+1), D)
        ∂ycoeffs = ntuple(_ -> zeros(degree+1, degree+1), D)
        differentiate!(∂xcoeffs, ∂ycoeffs, DiffPlan{Float64}(degree), phasecoeffs)

        phasepoints = [Vector{Float64}(undef, getpaduanum(degree)) for _ in 1:D]
        iplan = InvPaduaTransformPlan{Float64}(degree)
        invpaduatransform!(phasepoints, iplan, phasecoeffs)

        Ω(v, t, p) = CanonicalSymplecticMatrix(D)

        plan = CallIntPlan{Float64}(D, degree)

        intcoeffs = zeros(degree+1, degree+1)

        getintegrand!(intcoeffs, plan, Ω, phasepoints, 0, nothing, ∂xcoeffs, ∂ycoeffs)

        testcoeffs = zeros(degree+1, degree+1)
        testcoeffs[1, 1] = 0.25

        @test testcoeffs ≈ intcoeffs atol=1eps()
    end
end

@safetestset "compute! with ChebyshevPlan (OOP)" begin
    using ..Chebyshev.PaduaTransforms
    using ..Chebyshev: ChebyshevPlan, compute!, DiffPlan
    using PoincareInvariants.SecondPoincareInvariants
    using PoincareInvariants.CanonicalSymplecticStructures


    D = 12
    N = 200
    Ω(v, t, p) = CanonicalSymplecticMatrix(D)
    pointnum = nextpaduanum(N)
    degree = getdegree(pointnum)
    plan = ChebyshevPlan{Float64}(Ω, D, pointnum)

    testphasecoeffs = [zeros(degree+1, degree+1) for _ in 1:D]
    test∂x = [zeros(degree+1, degree+1) for _ in 1:D]
    test∂y = [zeros(degree+1, degree+1) for _ in 1:D]
    testintcoeffs = zeros(degree+1, degree+1)

    phasepoints = getpaduapoints(degree) do x, y
        ntuple(D) do i
            i == 1 && return x
            i == D ÷ 2 + 1 && return y
            return 0
        end
    end

    pinv = SecondPoincareInvariant{Float64}(Ω, D, pointnum, plan)
    @test compute!(pinv, phasepoints, 0, nothing) ≈ 4 atol=20eps()

    testphasecoeffs[1][1, 2] = 1  # 1 x
    testphasecoeffs[D ÷ 2 + 1][2, 1] = 1 # 1 y

    test∂x[1][1, 1] = 1
    test∂y[D ÷ 2 + 1][1, 1] = 1


    for d in 1:D
        @test maximum(abs, plan.phasecoeffs[d] .- testphasecoeffs[d]) / eps() < 10
        @test maximum(abs, plan.∂x[d] .- test∂x[d]) / eps() < 50
        @test maximum(abs, plan.∂y[d] .- test∂y[d]) / eps() < 50
    end

    testintcoeffs[1, 1] = 1
    @test maximum(abs, plan.intcoeffs .- testintcoeffs) / eps() < 50
end

@safetestset "getpoints, getpointspec and getpointnum" begin
    using ..Chebyshev.PaduaTransforms
    using ..Chebyshev: ChebyshevPlan, getpoints, getpointnum
    using PoincareInvariants.SecondPoincareInvariants: getpointspec

    for N in [10, 321, 2178], T in [Float32, Float64]
        @test getpointspec(N, ChebyshevPlan) == nextpaduanum(N)
        @test getpointspec((N, 2N), ChebyshevPlan) == nextpaduanum(N * 2N)

        @test getpointnum(N, ChebyshevPlan) == nextpaduanum(N)
        @test getpointnum((7, N), ChebyshevPlan) == nextpaduanum(7 * N)

        f(x, y) = 5 * x, cos(x*y), x+y
        rescale(x, y) = ((x, y) .+ 1) ./ 2

        pnts = getpoints(f, T, N, ChebyshevPlan)
        testpnts = getpaduapoints(T, nextdegree(N)) do x, y
            px, py = rescale(x, y)
            return f(px, py)
        end

        for i in 1:3
            @test pnts[i] ≈ testpnts[i]
        end
    end
end
