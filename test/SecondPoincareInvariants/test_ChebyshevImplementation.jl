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

    getintegrand!(intcoeffs, plan, Ω, phasepoints, 0, nothing, phasecoeffs, ∂xcoeffs, ∂ycoeffs)

    @test intcoeffs ≈ [-72   -82  -30;
                       -82  -432    0;
                       -30     0    0] atol=5eps()
end
