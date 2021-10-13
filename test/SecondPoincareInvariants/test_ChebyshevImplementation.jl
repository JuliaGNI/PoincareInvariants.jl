@safetestset "PaduaTransforms" begin include("test_PaduaTransforms.jl") end

@safetestset "Differentiation" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation:
        DiffPlan, differentiate!

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

    coeffsarr = rand(4, 4, 6)
    ∂xarr = zeros(4, 4, 6)
    ∂yarr = zeros(4, 4, 6)

    differentiate!(∂xarr, ∂yarr, P, coeffsarr)

    for i in 1:6
        differentiate!(∂x, ∂y, P, coeffsarr[:, :, i])
        @test ∂xarr[:, :, i] == ∂x
        @test ∂yarr[:, :, i] == ∂y
    end
end

@safetestset "Integration" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation:
        getintegrator, integrate

    # integrating odd polynomials over symmetric boundary conditions gives 0
    @test all(getintegrator(100)[2:2:end] .== 0)
    @test all(getintegrator(999)[2:2:end] .== 0)

    coeffs = [1 4 7 0;
              2 5 8 0;
              3 6 9 0;
		      0 0 0 0]

    @test integrate(coeffs, getintegrator(3)) ≈ 1 * 4 + 3 * -4/3 + 7 * -4/3 + 9 * 4/9 atol=5eps()
end

@safetestset "getintegrand! with OOPIntPlan" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.PaduaTransforms
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation:
        DiffPlan, differentiate!, OOPIntPlan, getintegrand!
    using PoincareInvariants.CanonicalSymplecticStructures

    @testset "6 Points in 2 Dimensions" begin
        degree = 2
        D = 2
        phasecoeffs = cat(
            Float64[1 2 3;
                    4 5 0;
                    6 0 0],
            Float64[1 4 6;
                    2 5 0;
                    3 0 0]
        ; dims=3)

        ∂xcoeffs, ∂ycoeffs = differentiate!(similar(phasecoeffs), similar(phasecoeffs),
            DiffPlan{Float64}(degree), phasecoeffs)

        phasepoints = invpaduatransform!(zeros(getpaduanum(degree), D),
            InvPaduaTransformPlan{Float64}(degree), phasecoeffs)

        Ω(v, t, p) = [0 -1; 1  0]

        plan = OOPIntPlan{Float64}(D, degree)

        intcoeffs = zeros(degree+1, degree+1)

        getintegrand!(intcoeffs, plan, Ω, phasepoints, 0, nothing, ∂xcoeffs, ∂ycoeffs)

        @test intcoeffs ≈ [-72   -82  -30;
                           -82  -432    0;
                           -30     0    0] atol=5eps()
    end

    @testset "$N Points in $D Dimensions" for D in [4, 20], N in [53, 11223]
        pointnum = nextpaduanum(N)
        degree = getdegree(pointnum)

        phasecoeffs = zeros(degree+1, degree+1, D)
        phasecoeffs[1, 1, 1] = 0.5  # 0.5
        phasecoeffs[1, 2, 1] = 0.5  # 0.5 x
        phasecoeffs[1, 1, D ÷ 2 + 1] = 0.5  # 0.5
        phasecoeffs[2, 1, D ÷ 2 + 1] = 0.5  # 0.5 y

        ∂xcoeffs, ∂ycoeffs = differentiate!(similar(phasecoeffs), similar(phasecoeffs),
            DiffPlan{Float64}(degree), phasecoeffs)

        phasepoints = invpaduatransform!(zeros(getpaduanum(degree), D),
            InvPaduaTransformPlan{Float64}(degree), phasecoeffs)

        Ω(v, t, p) = CanonicalSymplecticMatrix(D)

        plan = OOPIntPlan{Float64}(D, degree)

        intcoeffs = zeros(degree+1, degree+1)

        getintegrand!(intcoeffs, plan, Ω, phasepoints, 0, nothing, ∂xcoeffs, ∂ycoeffs)

        testcoeffs = zeros(degree+1, degree+1)
        testcoeffs[1, 1] = 0.25

        @test intcoeffs ≈ intcoeffs atol=1eps()
    end
end

@safetestset "compute! with ChebyshevPlan (OOP)" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.PaduaTransforms
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation:
        ChebyshevPlan, compute!, DiffPlan
    using PoincareInvariants.CanonicalSymplecticStructures


    @testset "$N Points in $D dimensions" for D in [4, 12], N in [20, 1234]
        Ω(v, t, p) = CanonicalSymplecticMatrix(D)
        plan = ChebyshevPlan{Float64}(Ω, D, N, Val(false))

        pointnum = nextpaduanum(N)
        degree = getdegree(pointnum)

        @test plan.degree == degree

        @test plan.paduaplan isa PaduaTransformPlan
        @test plan.paduaplan.degree == degree

        testphasecoeffs = zeros(degree+1, degree+1, D)
        @test plan.phasecoeffs == testphasecoeffs

        @test plan.diffplan isa DiffPlan
        @test size(plan.diffplan.D) == (degree+1, degree+1)

        test∂x = zeros(degree+1, degree+1, D)
        @test size(plan.∂x) == size(test∂x)

        test∂y = zeros(degree+1, degree+1, D)
        @test size(plan.∂y) == size(test∂y)

        testintcoeffs = zeros(degree+1, degree+1)
        @test size(plan.intcoeffs) == size(testintcoeffs)
        @test size(plan.integrator) == (degree+1,)

        phasepoints = getpaduapoints(degree) do x, y
            ntuple(D) do i
                i == 1 && return x
                i == D ÷ 2 + 1 && return y
                return 0
            end
        end

        @test compute!(plan, Ω, phasepoints, 0, nothing) ≈ 4 atol=10eps()

        testphasecoeffs[1, 2, 1] = 1  # 1 x
        testphasecoeffs[2, 1, D ÷ 2 + 1] = 1 # 1 y
        @test maximum(abs, plan.phasecoeffs .- testphasecoeffs) / eps() < 10

        test∂x[1, 1, 1] = 1
        @test maximum(abs, plan.∂x .- test∂x) / eps() < 50

        test∂y[1, 1, D ÷ 2 + 1] = 1
        @test maximum(abs, plan.∂y .- test∂y) / eps() < 50

        testintcoeffs[1, 1] = 1
        @test maximum(abs, plan.intcoeffs .- testintcoeffs) / eps() < 50
    end
end

@safetestset "getpoints and getpointnum" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.PaduaTransforms
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation:
        ChebyshevPlan, getpoints, getpointnum

    Ω(v, t, p) = [0 -1; 1 0]
    plan = ChebyshevPlan{Float64}(Ω, 2, 11, Val(false))

    @test getpointnum(plan) == 15
    @test getpoints(plan) == map(getpaduapoints(4)) do v
        (v .+ 1) ./ 2
    end
end
