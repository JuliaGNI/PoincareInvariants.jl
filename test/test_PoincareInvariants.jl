@safetestset "Properties and Basic Interface" begin
    @safetestset "FirstPoincareInvariant" begin
        using PoincareInvariants

        D = 6; N = 123; θ = canonical_one_form; P = FirstFinDiffPlan
        pinv = FirstPoincareInvariant{Float64, D}(θ, N, P)

        @test getdim(pinv) == D
        @test getplan(pinv) isa P
        @test N <= getpointnum(pinv) <= 2N
        @test getpointspec(pinv) == getpointspec(N, typeof(getplan(pinv)))
        @test getform(pinv) === θ

        pnts = getpoints(pinv) do x
            x, x, x, 0, 1, 2
        end
        @test pnts isa Matrix{Float64}
        @test size(pnts) == (getpointnum(pinv), D)

        @test size(getpoints(pinv)::Vector{Float64}) == (getpointnum(pinv),)
    end

    @safetestset "SecondPoincareInvariant" begin
        using PoincareInvariants

        D = 4; N = 321; ω = canonical_two_form; P = SecondFinDiffPlan
        pinv = SecondPoincareInvariant{Float64, D}(ω, N, P)

        @test getdim(pinv) == D
        @test getplan(pinv) isa P
        @test N <= getpointnum(pinv) <= 2N
        @test getpointspec(pinv) == getpointspec(N, typeof(getplan(pinv)))
        @test getform(pinv) === ω

        pnts = getpoints(pinv) do x, y
            x, x, y, y
        end
        @test pnts isa Matrix{Float64}
        @test size(pnts) == (getpointnum(pinv), D)

        @test size(getpoints(pinv)::Matrix{Float64}) == (getpointnum(pinv), 2)
    end
end

@safetestset "Basic Usage" begin
    @safetestset "FirstPoincareInvariant" begin
        using PoincareInvariants
        using DoubleFloats

        D = 2

        Ts = [
            Float64, Float64, Double64, Float64, Double64,
            Float64, Double64, Double64,
            Double64, Float64, Float64
        ]

        Ns = [
            500, 501, 518, 519, 520,
            541, 542, 543,
            566, 567, 568
        ]

        first_pinvs = [
            FirstPoincareInvariant{Ts[1], D}(canonical_one_form, Ns[1]),
            FirstPoincareInvariant{Ts[2], D}(canonical_one_form, Ns[2], FirstFinDiffPlan),
            FirstPoincareInvariant{Ts[3], D}(canonical_one_form, Ns[3], FirstFinDiffPlan),
            FirstPoincareInvariant{Ts[4], D}(canonical_one_form, Ns[4], FirstFourierPlan),
            FirstPoincareInvariant{Ts[5], D}(canonical_one_form, Ns[5], FirstFourierPlan),

            FirstPI{Ts[6], D}(canonical_one_form, Ns[6]),
            FirstPI{Ts[7], D}(canonical_one_form, Ns[7], FirstFinDiffPlan),
            FirstPI{Ts[8], D}(canonical_one_form, Ns[8], FirstFourierPlan),

            CanonicalFirstPI{Ts[9], D}(Ns[9]),
            CanonicalFirstPI{Ts[10], D}(Ns[10], FirstFinDiffPlan),
            CanonicalFirstPI{Ts[11], D}(Ns[11], FirstFourierPlan),
        ]

        to_circle_line(θ) = (cospi(-2θ), sinpi(-2θ))

        @testset "FirstPoincareInvariant $i" for (i, pinv) in enumerate(first_pinvs)
            T = Ts[i]
            N = Ns[i]

            pnts = getpoints(to_circle_line, pinv)
            @test pnts isa Matrix{T}
            @test size(pnts, 2) == D
            @test size(pnts, 1) == getpointnum(pinv)
            @test N ≤ getpointnum(pinv) ≤ 2N

            I = compute!(pinv, pnts, 0.1, nothing)
            I == compute!(pinv, pnts, 0.1)
            @test I isa T
            @test π ≈ I atol=1e-4
        end
    end

    @safetestset "SecondPoincareInvariant" begin
        using PoincareInvariants
        using DoubleFloats

        D = 2

        Ts = [
            Float64, Float64, Double64, Float64, Double64,
            Float64, Float64, Double64,
            Float64, Double64, Float64,
            Double64, Float64, Float64
        ]

        ps = [
            100, 101, 118, (11, 12), 140,
            162, 173, 196,
            222, 230, 244,
            253, 289, (19, 16)
        ]

        second_pinvs = [
            SecondPoincareInvariant{Ts[1], D}(canonical_two_form, ps[1]),
            SecondPoincareInvariant{Ts[2], D}(canonical_two_form, ps[2], SecondChebyshevPlan),
            SecondPoincareInvariant{Ts[3], D}(canonical_two_form, ps[3], SecondChebyshevPlan),
            SecondPoincareInvariant{Ts[4], D}(canonical_two_form, ps[4], SecondFinDiffPlan),
            SecondPoincareInvariant{Ts[5], D}(canonical_two_form, ps[5], SecondFinDiffPlan),

            SecondPI{Ts[6], D}(canonical_two_form, ps[6]),
            SecondPI{Ts[7], D}(canonical_two_form, ps[7], SecondChebyshevPlan),
            SecondPI{Ts[8], D}(canonical_two_form, ps[8], SecondFinDiffPlan),

            SecondPI{Ts[9], D}(CanonicalSymplecticMatrix{Float64}(D), ps[9]),
            SecondPI{Ts[10], D}(CanonicalSymplecticMatrix{Double64}(D), ps[10], SecondChebyshevPlan),
            SecondPI{Ts[11], D}(CanonicalSymplecticMatrix{Double64}(D), ps[11], SecondFinDiffPlan),

            CanonicalSecondPI{Ts[12], D}(ps[12]),
            CanonicalSecondPI{Ts[13], D}(ps[13], SecondChebyshevPlan),
            CanonicalSecondPI{Ts[14], D}(ps[14], SecondFinDiffPlan)
        ]

        to_wavy_area(x, y) = (x, y + sinpi(2x))

        @testset "SecondPoincareInvariant $i" for (i, pinv) in enumerate(second_pinvs)
            T = Ts[i]
            psi = ps[i]
            N = getpointnum(psi)

            pnts = getpoints(to_wavy_area, pinv)
            @test pnts isa Matrix{T}
            @test size(pnts, 2) == D
            @test size(pnts, 1) == getpointnum(pinv)
            @test N ≤ getpointnum(pinv) ≤ 2N

            I = compute!(pinv, pnts, 0.1, nothing)
            I == compute!(pinv, pnts, 0.1)
            @test I isa T
            @test 1 ≈ I atol=50eps()
        end
    end
end

@safetestset "Consistency Between Implementations" begin
    @safetestset "FirstPoincareInvariant" begin
        using PoincareInvariants

        D = 6
        θ(z, t, p) = [-z[4], -z[5], -p * z[6], t * z[1], p * z[2], z[3]]

        function f(ϕ)
            y, x = sincospi(2ϕ)
            return 0.5 * x, x, 3 * x, 0.1 * y, 0.85 * y, 2 * y
        end

        Idefault = let pinv = FirstPI{Float64, D}(θ, 1_000)
            compute!(pinv, getpoints(f, pinv), 0.3, 0.1)
        end

        Idiff = let pinv = FirstPI{Float64, D}(θ, 1_000, FirstFinDiffPlan)
            compute!(pinv, getpoints(f, pinv), 0.3, 0.1)
        end

        Ifreq = let pinv = FirstPI{Float64, D}(θ, 1_000, FirstFourierPlan)
            compute!(pinv, getpoints(f, pinv), 0.3, 0.1)
        end

        @test Ifreq ≈ Idefault atol=10eps()
        @test Ifreq ≈ Idiff atol=5e-4
    end

    @safetestset "SecondPoincareInvariant" begin
        using PoincareInvariants

        D = 6
        ω(z, t, p) = [
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

        Idefault = let pinv = SecondPI{Float64, D}(ω, 1_000)
            compute!(pinv, getpoints(f, pinv), 5, 11)
        end

        Icheb = let pinv = SecondPI{Float64, D}(ω, 1_000, SecondChebyshevPlan)
            compute!(pinv, getpoints(f, pinv), 5, 11)
        end

        Ifindiff = let pinv = SecondPI{Float64, D}(ω, (31, 33), SecondFinDiffPlan)
            compute!(pinv, getpoints(f, pinv), 5, 11)
        end

        @test Icheb ≈ Idefault atol=10eps()
        @test Icheb ≈ Ifindiff atol=5e-3
    end
end

@safetestset "Iterable of Points" begin
    using PoincareInvariants

    D = 4
    N = 398
    Nsteps = 13
    times = [2.0 + n * 0.15 for n in 0:Nsteps-1]

    f1(z, t, p) = canonical_one_form(z, t, p) .* t
    pi1 = FirstPI{Float64, D}(f1, N)

    f2(z, t, p) = canonical_two_form(z, t, p) .* t
    pi2 = SecondPI{Float64, D}(f2, N)

    p1(θ) = (cospi(2θ), cospi(2θ), sinpi(2θ), sinpi(2θ))
    p2(x, y) = (x, x, y, y)

    data1 = map(eachrow(getpoints(p1, pi1))) do row
        fill(ntuple(i -> row[i], Val(D)), length(times))
    end
    data2 = map(eachrow(getpoints(p2, pi2))) do row
        fill(ntuple(i -> row[i], Val(D)), length(times))
    end

    @test [-2π * t for t in times] ≈ compute!(pi1, data1, times) atol=1e-13
    @test [2 * t for t in times] ≈ compute!(pi2, data2, times, nothing) atol=1e-13
end

@safetestset "Linear Symplectic Maps" begin
    @safetestset "In 2D" begin
        using PoincareInvariants

        diffpi1 = CanonicalFirstPI{Float64, 2}(1_000, FirstFinDiffPlan)
        freqpi1 = CanonicalFirstPI{Float64, 2}(1_000, FirstFourierPlan)

        chebpi2 = CanonicalSecondPI{Float64, 2}(1_000, SecondChebyshevPlan)
        diffpi2 = CanonicalSecondPI{Float64, 2}(1_000, SecondFinDiffPlan)

        A = [1 5;
             0 1]
        B = [2   0;
             0 0.5]

        # circle (negative orientation)
        f1(θ) = A * B * [cospi(2θ), sinpi(2θ)]

        # half circle with radius one
        f2(r, θ) = A * B * [r*cospi(θ), r*sinpi(θ)]

        @test -π ≈ compute!(diffpi1, getpoints(f1, diffpi1)) atol=5e-5
        @test -π ≈ compute!(freqpi1, getpoints(f1, freqpi1)) atol=10eps()

        @test π/2 ≈ compute!(chebpi2, getpoints(f2, chebpi2)) atol=10eps()
        @test π/2 ≈ compute!(diffpi2, getpoints(f2, diffpi2)) atol=5e-3
    end

    @safetestset "In 8D" begin
        using PoincareInvariants

        diffpi1 = CanonicalFirstPI{Float64, 8}(1_000, FirstFinDiffPlan)
        freqpi1 = CanonicalFirstPI{Float64, 8}(1_000, FirstFourierPlan)

        chebpi2 = CanonicalSecondPI{Float64, 8}(1_000, SecondChebyshevPlan)
        diffpi2 = CanonicalSecondPI{Float64, 8}(1_000, SecondFinDiffPlan)

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
               0   0    0    0   0 0.2    0 -0.4;
               0   0    0    0 0.4   0 -0.2    0]

        # test if symplectic
        @test transpose(A * B * C) * (-A) * (A * B * C) ≈ -A

        # half circle with positive orientation
        f1(θ) = A * B * C * [0, cospi(-θ), 0, 0, 0, sinpi(-θ), 0, 0]

        # half circle with radius one and negative orientation
        f2(r, θ) = A * B * C * [0, r*cospi(θ), 0, 0, 0, -r*sinpi(θ), 0, 0]

        @test π/2 ≈ compute!(diffpi1, getpoints(f1, diffpi1)) atol=1e-2
        @test π/2 ≈ compute!(freqpi1, getpoints(f1, freqpi1)) atol=5e-3

        @test -π/2 ≈ compute!(chebpi2, getpoints(f2, chebpi2)) atol=10eps()
        @test -π/2 ≈ compute!(diffpi2, getpoints(f2, diffpi2)) atol=5e-3
    end
end

@safetestset "Henon Map" begin
    using PoincareInvariants

    init1(θ) = (cospi(2θ), cospi(2θ), 0, -sinpi(2θ))
    I1 = π

    init2(x, y) = 2 .* (x - 0.5, -x + 0.5, y - 0.5, -y + 0.5)
    I2 = 2^2 * 2

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

    diffpi1 = CanonicalFirstPI{Float64, 4}(1_000, FirstFinDiffPlan)
    freqpi1 = CanonicalFirstPI{Float64, 4}(1_000, FirstFourierPlan)

    chebpi2 = CanonicalSecondPI{Float64, 4}(1_000, SecondChebyshevPlan)
    diffpi2 = CanonicalSecondPI{Float64, 4}(1_000, SecondFinDiffPlan)

    f1(x) = init1(x) |> f
    @test I1 ≈ compute!(diffpi1, getpoints(f1, diffpi1)) atol=5e-5
    @test I1 ≈ compute!(freqpi1, getpoints(f1, freqpi1)) atol=5e-15

    g1(x, y) = init2(x, y) |> f
    @test I2 ≈ compute!(chebpi2, getpoints(g1, chebpi2)) atol=1e-14
    @test I2 ≈ compute!(diffpi2, getpoints(g1, diffpi2)) atol=5e-13

    f2(x) = init1(x) |> f |> f
    @test I1 ≈ compute!(diffpi1, getpoints(f2, diffpi1)) atol=5e-5
    @test I1 ≈ compute!(freqpi1, getpoints(f2, freqpi1)) atol=5e-13

    g2(x, y) = init2(x, y) |> f |> f
    @test I2 ≈ compute!(chebpi2, getpoints(g2, chebpi2)) atol=1e-11
    @test I2 ≈ compute!(diffpi2, getpoints(g2, diffpi2)) atol=5e-12

    f3(x) = init1(x) |> f |> f |> f
    @test I1 ≈ compute!(freqpi1, getpoints(f3, freqpi1)) atol=1e-7

    g3(x, y) = init2(x, y) |> f |> f |> f
    @test I2 ≈ compute!(chebpi2, getpoints(g3, chebpi2)) atol=5e-6
end

@safetestset "Wavy Map" begin
    using PoincareInvariants

    function wavy(x, y)
        x, y = x, y + 2 * sinpi(x)
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

    init1(θ) = (cospi(2θ), sinpi(2θ), cospi(2θ), sinpi(2θ), 1, sinpi(2θ))
    I1 = -2π

    init2(x, y) = (2x, 2x-1, 2x-2, 2y-1, 2y-1, 2y-1)
    I2 = 3 * 2^2

    diffpi1 = CanonicalFirstPI{Float64, 6}(1_000, FirstFinDiffPlan)
    freqpi1 = CanonicalFirstPI{Float64, 6}(1_000, FirstFourierPlan)

    chebpi2 = CanonicalSecondPI{Float64, 6}(1_000, SecondChebyshevPlan)
    diffpi2 = CanonicalSecondPI{Float64, 6}(1_000, SecondFinDiffPlan)

    f1(x) = init1(x) |> f
    @test I1 ≈ compute!(diffpi1, getpoints(f1, diffpi1)) atol=5e-5
    @test I1 ≈ compute!(freqpi1, getpoints(f1, freqpi1)) atol=1e-15

    g1(x, y) = init2(x, y) |> f
    @test I2 ≈ compute!(chebpi2, getpoints(g1, chebpi2)) atol=1e-14
    @test I2 ≈ compute!(diffpi2, getpoints(g1, diffpi2)) atol=5e-15

    f2(x) = init1(x) |> f |> f
    @test I1 ≈ compute!(diffpi1, getpoints(f2, diffpi1)) atol=1e-3
    @test I1 ≈ compute!(freqpi1, getpoints(f2, freqpi1)) atol=5e-14

    g2(x, y) = init2(x, y) |> f |> f
    @test I2 ≈ compute!(chebpi2, getpoints(g2, chebpi2)) atol=5e-3
    @test I2 ≈ compute!(diffpi2, getpoints(g2, diffpi2)) atol=1e-3

    f3(x) = init1(x) |> f |> f |> f
    @test I1 ≈ compute!(diffpi1, getpoints(f3, diffpi1)) atol=5e-2
    @test I1 ≈ compute!(freqpi1, getpoints(f3, freqpi1)) atol=1e-13

    g3(x, y) = init2(x, y) |> f |> f |> f
    # neither of the two methods manage this
end

@safetestset "Pendulum" begin
    using OrdinaryDiffEq
    using RecursiveArrayTools: ArrayPartition
    using PoincareInvariants

    # ArrayPartition is a pain to work with!
    prob = SecondOrderODEProblem((p, θ, params, t) -> [-sin(θ[1])], 0.0, 0.0, (0.0, 2.0))
    dt = 0.1
    pf(prob, i, repeat) = remake(prob; u0 = ArrayPartition((prob.u0[1:1], prob.u0[2:2])))

    I1 = 3π
    pi1 = CanonicalFirstPI{Float64, 2}(1_000)
    prob1 = PIEnsembleProblem(ϕ -> (sinpi(2ϕ), 3 * cospi(2ϕ)), prob, pi1; prob_func=pf)
    sol1 = solve(prob1, SymplecticEuler(), EnsembleSerial(); adaptive=false, dt=dt)

    I2 = 16
    pi2 = CanonicalSecondPI{Float64, 2}(1_000)
    prob2 = PIEnsembleProblem((x, y) -> 4 .* (x, y) .- 2, prob, pi2; prob_func=pf)
    sol2 = solve(prob2, SymplecticEuler(), EnsembleSerial(); adaptive=false, dt=dt)

    @test maximum(abs, I1 .- compute!(pi1, sol1)) < 1e-14
    @test maximum(abs, I2 .- compute!(pi2, sol2)) < 1e-10
end
