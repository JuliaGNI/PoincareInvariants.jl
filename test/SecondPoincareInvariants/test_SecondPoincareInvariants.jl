@safetestset "Chebyshev Implementation" begin include("test_Chebyshev.jl") end

@safetestset "FiniteDifferences Implementation" begin
    include("test_FiniteDifferences.jl")
end

@safetestset "One by One Square" begin
    using PoincareInvariants
    using PoincareInvariants.SecondPoincareInvariants.Chebyshev: ChebyshevPlan
    using PoincareInvariants.SecondPoincareInvariants.FiniteDifferences: FiniteDiffPlan

    D = 2
    Ω(z, t, p) = CanonicalSymplecticMatrix(D)

    let pinv = SecondPoincareInvariant{Float64}(Ω, D, 432)
        I = compute!(pinv, getpoints(pinv), 0, nothing)
        @test abs(1 - I) / eps() < 10
    end

    let pinv = SecondPoincareInvariant{Float64}(Ω, D, 567, ChebyshevPlan)
        I = compute!(pinv, getpoints(pinv), 0, nothing)
        @test abs(1 - I) / eps() < 10
    end

    let pinv = SecondPoincareInvariant{Float64}(Ω, D, (35, 75), FiniteDiffPlan)
        I = compute!(pinv, getpoints(pinv), 0, nothing)
        @test abs(1 - I) / eps() < 10
    end
end

@safetestset "Consistency Between Implementations" begin
    using PoincareInvariants
    using PoincareInvariants.SecondPoincareInvariants.Chebyshev: ChebyshevPlan
    using PoincareInvariants.SecondPoincareInvariants.FiniteDifferences: FiniteDiffPlan

    D = 6
    Ω(z, t, p) = [
            0     0     0  z[1]  z[2]  z[3]
            0     0     0  z[4]  z[5]  z[6]
            0     0     0     0     0     0
        -z[1] -z[4]     0     0     0     0
        -z[2] -z[5]     0     0     0     0
        -z[3] -z[6]     0     0     0     0
    ]

    f(x, y) = (
        exp(x) * y,
        cospi(x * y),
        x^2 - y^3 + x*y,
        5.0,
        exp(-(x^2 + y^2)),
        x + y
    )

    Idefault = let pinv = PI2{Float64}(Ω, D, 10_000)
        compute!(pinv, getpoints(f, pinv), 0, nothing)
    end

    Icheb = let pinv = PI2{Float64}(Ω, D, 10_000, ChebyshevPlan)
        compute!(pinv, getpoints(f, pinv), 0, nothing)
    end

    Ifindiff = let pinv = PI2{Float64}(Ω, D, (100, 100), FiniteDiffPlan)
        compute!(pinv, getpoints(f, pinv), 0, nothing)
    end

    @test Icheb ≈ Idefault rtol=10eps()
    @test Icheb ≈ Ifindiff rtol=1e-3
end

@safetestset "Free Particle" begin
    using PoincareInvariants

    function free_particle!(points, t)
        mid = length(points) ÷ 2
        for i in 1:mid
            points[i] .+= points[mid+i] .* t
        end
    end

    @testset "Square" begin
        D = 8
        N = 1_000
        Ω(z, t, p) = CanonicalSymplecticMatrix(D)
        pinv = SecondPoincareInvariant{Float64}(Ω, D, N)

        phasepoints = getpoints(pinv) do x, y
            (x, 0, 0, 0, y, 0, 0, 0)
        end

        @test abs(1 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 10

        free_particle!(phasepoints, 10)
        @test abs(1 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 20

        free_particle!(phasepoints, 100)
        @test abs(1 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 100
    end

    @testset "Quarter Circle" begin
        D = 4
        N = 2_000
        Ω(z, t, p) = CanonicalSymplecticMatrix(D)
        pinv = SecondPoincareInvariant{Float64}(Ω, D, N)

        phasepoints = getpoints(pinv) do r, θ
            (0, r * cos(θ * π/2), 0, r * sin(θ * π/2))
        end

        @test abs(π/4 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 10

        free_particle!(phasepoints, 10)
        @test abs(π/4 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 20

        free_particle!(phasepoints, 100)
        @test abs(π/4 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 200
    end

    @testset "Half Sphere in 4D" begin
        D = 4
        N = 1500
        Ω(z, t, p) = CanonicalSymplecticMatrix(D)
        pinv = SecondPoincareInvariant{Float64}(Ω, D, N)

        phasepoints = getpoints(pinv) do θ, ϕ
            sinθ, cosθ = sincospi(θ)
            sinϕ, cosϕ = sincospi(ϕ)
            return (sinθ * cosϕ, sinθ * sinϕ, cosθ, -1)
        end

        # invariant is -π
        # TODO: do the calculation by hand and double check it

        @test abs(-π - compute!(pinv, phasepoints, 0, nothing)) / eps() < 15

        free_particle!(phasepoints, 10)
        @test abs(-π - compute!(pinv, phasepoints, 0, nothing)) / eps() < 15

        free_particle!(phasepoints, 100)
        @test abs(-π - compute!(pinv, phasepoints, 0, nothing)) / eps() < 150
    end
end
