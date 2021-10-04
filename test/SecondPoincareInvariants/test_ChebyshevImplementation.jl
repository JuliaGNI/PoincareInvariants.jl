@safetestset "Padua" begin include("test_Padua.jl") end

@safetestset "PaduaSetup and paduatransform!" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation.Padua
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation:
        PaduaSetup, paduatransform!

    degree = 20
    dims = 5
    plan = PaduaTransformPlan{Float64}(degree)

    vals = rand(getpaduanum(degree), dims)
    out = Matrix{Float64}(undef, getpaduanum(degree), dims)

    paduatransform!(out, plan, vals, Val(false))

    setup = PaduaSetup{Float64}(dims, degree)
    coeffs = paduatransform!(setup, vals, Val(false))

    @test coeffs == out
end

@safetestset "DiffSetup and differentiate!" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation:
        DiffSetup, differentiate!

    setup = DiffSetup{Float64}(3, 2)

    @test setup.Dx == [0 0 1 0  0  0; # U00
                       0 0 0 0 0.5 0; # U01
                       0 0 0 0  0  2] # U10

    @test setup.Dy == [0 1 0 0  0  0; # U00
                       0 0 0 2  0  0; # U01
                       0 0 0 0 0.5 0] # U10

    coeffs = rand(6, 3)

    ∂x, ∂y = differentiate!(setup, coeffs)

    @test setup.∂xcoeffs === ∂x == setup.Dx * coeffs
    @test size(setup.∂xcoeffs) == (3, 3)

    @test setup.∂ycoeffs === ∂y == setup.Dy * coeffs
    @test size(setup.∂ycoeffs) == (3, 3)
end

@safetestset "C12Setup and C12convert!" begin
    using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation:
        C12Setup, C12convert!

    setup = C12Setup{Float64}(5, 2)

    @test setup.C12 == [ 1  0   0  -0.5  0  -0.5; # U00
                         0 0.5  0    0   0    0 ; # U01
                         0  0  0.5   0   0    0 ] # U10

    coeffs = rand(6, 5)

    out = C12convert!(setup, coeffs)

    @test setup.coeffs === out == setup.C12 * coeffs
    @test size(setup.coeffs) == (3, 5)
end
