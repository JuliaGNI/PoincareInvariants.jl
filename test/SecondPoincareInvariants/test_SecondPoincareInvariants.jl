@safetestset "Chebyshev Implementation" begin include("test_Chebyshev.jl") end

@safetestset "FiniteDifferences Implementation" begin
    include("test_FiniteDifferences.jl")
end

@safetestset "Interface" begin
    @safetestset "Basic properties" begin
        using PoincareInvariants

        D = 6; N = 123
        Ω(z, t, p) = CanonicalSymplecticMatrix(D)
        pinv = SecondPoincareInvariant{Float64, D}(Ω, N)

        @test getdim(pinv) == 6
        @test N <= getpointnum(pinv) <= 2N
        @test getpointspec(pinv) == getpointspec(N, typeof(pinv.plan)) == pinv.pointspec
        @test getform(pinv) === pinv.Ω === Ω
    end
end

@safetestset "One by One Square" begin
    using PoincareInvariants
    using PoincareInvariants.SecondPoincareInvariants.Chebyshev: ChebyshevPlan
    using PoincareInvariants.SecondPoincareInvariants.FiniteDifferences: FiniteDiffPlan

    D = 2
    Ω(z, t, p) = CanonicalSymplecticMatrix(D)

    let pinv = SecondPoincareInvariant{Float64, D}(Ω, 432)
        I = compute!(pinv, getpoints(pinv), 0, nothing)
        @test abs(1 - I) / eps() < 10
    end

    let pinv = SecondPoincareInvariant{Float64, D}(Ω, 567, ChebyshevPlan)
        I = compute!(pinv, getpoints(pinv), 0, nothing)
        @test abs(1 - I) / eps() < 10
    end

    let pinv = SecondPoincareInvariant{Float64, D}(Ω, (35, 75), FiniteDiffPlan)
        I = compute!(pinv, getpoints(pinv), 0, nothing)
        @test abs(1 - I) / eps() < 10
    end
end

@safetestset "Consistency Between Implementations" begin
    using PoincareInvariants
    using PoincareInvariants.SecondPoincareInvariants: ChebyshevPlan
    using PoincareInvariants.SecondPoincareInvariants: FiniteDiffPlan

    D = 6
    Ω(z, t, p) = [
            0     0     0  z[1]  z[2]  z[3]
            0     0     0  z[4]  z[5]  z[6]
            0     0     0     t     p     1
        -z[1] -z[4]    -t     0     0     0
        -z[2] -z[5]    -p     0     0     0
        -z[3] -z[6]    -1     0     0     0
    ]

    f(x, y) = (
        exp(x) * y,
        cospi(x * y),
        x^2 - y^3 + x*y,
        5.0,
        exp(-(x^2 + y^2)),
        x + y
    )

    Idefault = let pinv = PI2{Float64, D}(Ω, 1_000)
        compute!(pinv, getpoints(f, pinv), 5, 11)
    end

    Icheb = let pinv = PI2{Float64, D}(Ω, 1_000, ChebyshevPlan)
        compute!(pinv, getpoints(f, pinv), 5, 11)
    end

    Ifindiff = let pinv = PI2{Float64, D}(Ω, (31, 33), FiniteDiffPlan)
        compute!(pinv, getpoints(f, pinv), 5, 11)
    end

    @test Icheb ≈ Idefault atol=10eps()
    @test Icheb ≈ Ifindiff atol=5e-3
end

@safetestset "Linear Symplectic Maps" begin
    @safetestset "In 2D" begin
        using PoincareInvariants
        using PoincareInvariants.SecondPoincareInvariants: FiniteDiffPlan

        Ω(z, t, p) = CanonicalSymplecticMatrix(2)
        pinv1 = PI2{Float64, 2}(Ω, 1_000)
        pinv2 = PI2{Float64, 2}(Ω, 1_000, FiniteDiffPlan)

        A = [1 5;
             0 1]
        B = [2   0;
             0 0.5]
        f(r, θ) = A * B * [r*cospi(θ), r*sinpi(θ)]
        # half circle with radius one

        @test abs(π/2 - compute!(pinv1, getpoints(f, pinv1), 0, nothing)) / eps() < 10
        @test abs(π/2 - compute!(pinv2, getpoints(f, pinv2), 0, nothing)) < 5e-3
    end

    @safetestset "In 8D" begin
        using PoincareInvariants
        using PoincareInvariants.SecondPoincareInvariants: FiniteDiffPlan

        Ω(z, t, p) = CanonicalSymplecticMatrix(8)
        pinv1 = PI2{Float64, 8}(Ω, 1_000)
        pinv2 = PI2{Float64, 8}(Ω, 1_000, FiniteDiffPlan)

        A = CanonicalSymplecticMatrix(8)

        B = [1 0 0 0 1 5 8 1;
             0 1 0 0 5 2 6 9;
             0 0 1 0 8 6 3 7;
             0 0 0 1 1 9 7 4;
             0 0 0 0 1 0 0 0;
             0 0 0 0 0 1 0 0;
             0 0 0 0 0 0 1 0;
             0 0 0 0 0 0 0 1]

        C = [1.0   0  2.0    0   0   0    0    0;
               0 2.0    0  1.0   0   0    0    0;
               0 1.0    0 -2.0   0   0    0    0;
             2.0   0 -1.0    0   0   0    0    0;
               0   0    0    0 0.2   0  0.4    0;
               0   0    0    0   0 0.4    0  0.2;
               0   0    0    0 0.2   0    0 -0.4;
               0   0    0    0 0.4   0 -0.2    0]

        f(r, θ) = A * B * C * [r*sinpi(θ), 0, 0, 0, r*cospi(θ), 0, 0, 0]
        # half circle with radius one and negative orientation

        @test abs(-π/2 - compute!(pinv1, getpoints(f, pinv1), 0, nothing)) / eps() < 10
        @test abs(-π/2 - compute!(pinv2, getpoints(f, pinv2), 0, nothing)) < 5e-3
    end
end

@safetestset "Henon Map" begin
    using PoincareInvariants
    using PoincareInvariants.SecondPoincareInvariants: FiniteDiffPlan

    Ω(z, t, p) = CanonicalSymplecticMatrix(4)
    pinv1 = PI2{Float64, 4}(Ω, 1_000)
    pinv2 = PI2{Float64, 4}(Ω, 1_000, FiniteDiffPlan)

    init(x, y) = (x, -x, y, -y)
    henon(q, p, a) = (p + a - q^2, -q)

    # henon plus symplectic intermixing of q1/p1 with q2/p2
    function f((q1, q2, p1, p2))
        q1, p1 = henon(q1, p1, 2)
        q2, p2 = henon(q2, p2, 3)
        return (q1 - p2, q2 - p1, p1, p2)
    end

    f1(x, y) = init(x, y) |> f
    @test abs(2 - compute!(pinv1, getpoints(f1, pinv1), 0, nothing)) / eps() < 10
    @test abs(2 - compute!(pinv2, getpoints(f1, pinv2), 0, nothing)) / eps() < 500

    f2(x, y) = init(x, y) |> f |> f
    @test abs(2 - compute!(pinv1, getpoints(f2, pinv1), 0, nothing)) / eps() < 50
    @test abs(2 - compute!(pinv2, getpoints(f2, pinv2), 0, nothing)) < 0.01

    f3(x, y) = init(x, y) |> f |> f |> f
    @test abs(2 - compute!(pinv1, getpoints(f3, pinv1), 0, nothing)) / eps() < 100
    @test abs(2 - compute!(pinv2, getpoints(f3, pinv2), 0, nothing)) < 0.5
end
