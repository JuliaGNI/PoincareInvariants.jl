@safetestset "Chebyshev Implementation" begin include("test_Chebyshev.jl") end

@safetestset "SecondPoincareInvariant OOP" begin
    using PoincareInvariants

    @testset "Points on Square" begin
        D = 10
        N = 500
        Ω(z, t, p) = CanonicalSymplecticMatrix(D)
        pinv = SecondPoincareInvariant{Float64}(Ω, D, N, Val(false))

        phasepoints = getpoints(pinv) do x, y
            ntuple(D) do i
                i == 1 && return x
                i == D ÷ 2 + 1 && return y
                return 0
            end
        end

        @test abs(1 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 10
    end
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
        pinv = SecondPoincareInvariant{Float64}(Ω, D, N, Val(false))

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
        pinv = SecondPoincareInvariant{Float64}(Ω, D, N, Val(false))

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
        pinv = SecondPoincareInvariant{Float64}(Ω, D, N, Val(false))

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
