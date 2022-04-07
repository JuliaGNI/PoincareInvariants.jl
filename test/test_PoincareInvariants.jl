@safetestset "Interface" begin
    @safetestset "Basic properties" begin
        using PoincareInvariants

        D = 6; N = 123
        Ω(z, t, p) = CanonicalSymplecticTwoForm(D)
        pinv = SecondPoincareInvariant{Float64, D}(Ω, N)

        @test getdim(pinv) == 6
        @test N <= getpointnum(pinv) <= 2N
        @test getpointspec(pinv) == getpointspec(N, typeof(pinv.plan)) == pinv.pointspec
        @test getform(pinv) === pinv.Ω === Ω
    end
end

@safetestset "One by One Square" begin
    using PoincareInvariants

    D = 2
    Ω(z, t, p) = CanonicalSymplecticTwoForm(D)

    let pinv = SecondPoincareInvariant{Float64, D}(Ω, 432)
        I = compute!(pinv, getpoints(pinv), 0, nothing)
        @test abs(1 - I) / eps() < 10
    end

    let pinv = SecondPoincareInvariant{Float64, D}(Ω, 567, SecondChebyshevPlan)
        I = compute!(pinv, getpoints(pinv), 0, nothing)
        @test abs(1 - I) / eps() < 10
    end

    let pinv = SecondPoincareInvariant{Float64, D}(Ω, (35, 75), SecondFinDiffPlan)
        I = compute!(pinv, getpoints(pinv), 0, nothing)
        @test abs(1 - I) / eps() < 10
    end
end

@safetestset "Consistency Between Implementations" begin
    using PoincareInvariants

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

    Idefault = let pinv = SecondPI{Float64, D}(Ω, 1_000)
        compute!(pinv, getpoints(f, pinv), 5, 11)
    end

    Icheb = let pinv = SecondPI{Float64, D}(Ω, 1_000, SecondChebyshevPlan)
        compute!(pinv, getpoints(f, pinv), 5, 11)
    end

    Ifindiff = let pinv = SecondPI{Float64, D}(Ω, (31, 33), SecondFinDiffPlan)
        compute!(pinv, getpoints(f, pinv), 5, 11)
    end

    @test Icheb ≈ Idefault atol=10eps()
    @test Icheb ≈ Ifindiff atol=5e-3
end

@safetestset "Linear Symplectic Maps" begin
    @safetestset "In 2D" begin
        using PoincareInvariants

        Ω(z, t, p) = CanonicalSymplecticTwoForm(2)
        chebpi = SecondPI{Float64, 2}(Ω, 1_000)
        diffpi = SecondPI{Float64, 2}(Ω, 1_000, SecondFinDiffPlan)

        A = [1 5;
             0 1]
        B = [2   0;
             0 0.5]
        f(r, θ) = A * B * [r*cospi(θ), r*sinpi(θ)]
        # half circle with radius one

        @test abs(π/2 - compute!(chebpi, getpoints(f, chebpi), 0, nothing)) / eps() < 10
        @test abs(π/2 - compute!(diffpi, getpoints(f, diffpi), 0, nothing)) < 5e-3
    end

    @safetestset "In 8D" begin
        using PoincareInvariants

        Ω(z, t, p) = CanonicalSymplecticTwoForm(8)
        chebpi = SecondPI{Float64, 8}(Ω, 1_000)
        diffpi = SecondPI{Float64, 8}(Ω, 1_000, SecondFinDiffPlan)

        A = CanonicalSymplecticTwoForm(8)

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

        @test abs(-π/2 - compute!(chebpi, getpoints(f, chebpi), 0, nothing)) / eps() < 10
        @test abs(-π/2 - compute!(diffpi, getpoints(f, diffpi), 0, nothing)) < 5e-3
    end
end

@safetestset "Henon Map" begin
    using PoincareInvariants

    init(x, y) = 2 .* (x - 0.5, -x + 0.5, y - 0.5, -y + 0.5)
    I = 2^2 * 2

    # Area preserving Henon map
    function henon(q, p)
        a = 1.4; b = 0.3
        q, p = q, 1 + p - a * q^2 # fold
        q, p = b * q, p / b # press
        q, p = -p, q  # rotate
        return q, p
    end

    # henon plus symplectic intermixing of q1/p1 with q2/p2
    function f((q1, q2, p1, p2))
        q1, p1 = henon(q1, p1)
        q2, p2 = henon(q2, p2)
        return (q1, q2, p1 - q2, p2 - q1)
    end

    Ω(z, t, p) = CanonicalSymplecticTwoForm(4)
    chebpi = SecondPI{Float64, 4}(Ω, 1_000)
    diffpi = SecondPI{Float64, 4}(Ω, 1_000, SecondFinDiffPlan)

    f1(x, y) = init(x, y) |> f
    @test abs(I - compute!(chebpi, getpoints(f1, chebpi), 0, nothing)) < 1e-14
    @test abs(I - compute!(diffpi, getpoints(f1, diffpi), 0, nothing)) < 5e-13

    f2(x, y) = init(x, y) |> f |> f
    @test abs(I - compute!(chebpi, getpoints(f2, chebpi), 0, nothing)) < 1e-11
    @test abs(I - compute!(diffpi, getpoints(f2, diffpi), 0, nothing)) < 5e-12

    f3(x, y) = init(x, y) |> f |> f |> f
    @test abs(I - compute!(chebpi, getpoints(f3, chebpi), 0, nothing)) < 5e-6
end

@safetestset "Wavy Map" begin
    using PoincareInvariants

    function wavy(x, y)
        x, y = x, y + sinpi(x)
        x, y = -y, x
        return x, y
    end

    # make wavy and mix dims
    function f((x1, x2, x3, y1, y2, y3))
        x1, y1 = wavy(x1, y1)
        x2, y2 = wavy(x2, y2)
        x3, y3 = wavy(x3, y3)
        return (x1 + y2, x2 + y1 + y3, x3 + y2, y1, y2, y3)
    end

    init(x, y) = (2x, 2x-1, 2x-2, 2y-1, 2y-1, 2y-1)
    I = 3 * 2^2

    chebpi = CanonicalSecondPI{Float64, 6}(1_000)
    diffpi = CanonicalSecondPI{Float64, 6}(1_000, SecondFinDiffPlan)

    f1(x, y) = init(x, y) |> f
    @test abs(I - compute!(chebpi, getpoints(f1, chebpi), 0, nothing)) < 1e-14
    @test abs(I - compute!(diffpi, getpoints(f1, diffpi), 0, nothing)) < 5e-15

    f2(x, y) = init(x, y) |> f |> f
    @test abs(I - compute!(chebpi, getpoints(f2, chebpi), 0, nothing)) < 5e-6
    @test abs(I - compute!(diffpi, getpoints(f2, diffpi), 0, nothing)) < 1e-4

    f3(x, y) = init(x, y) |> f |> f |> f
    @test abs(I - compute!(chebpi, getpoints(f3, chebpi), 0, nothing)) < 1
end
