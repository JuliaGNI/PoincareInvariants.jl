@safetestset "ChebyshevImplementation" begin include("test_ChebyshevImplementation.jl") end

@safetestset "SecondPoincareInvariant OOP" begin
    using PoincareInvariants

    @testset "$N Points on Square in $(D)D" for D in [2, 10], N in [10, 123, 4321]
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

    function free_particle!(state, δt)
        mid = length(state) ÷ 2
        state[1:mid] += state[mid+1:end] .* δt
    end

    @testset "$N Points on Square in $(D)D" for D in [2, 6], N in [500, 10_000]
        Ω(z, t, p) = CanonicalSymplecticMatrix(D)
        pinv = SecondPoincareInvariant{Float64}(Ω, D, N, Val(false))

        parampoints = getpoints(pinv)

        phasepoints = getpoints(pinv) do x, y
            ntuple(D) do i
                i == 1 && return x
                i == D ÷ 2 + 1 && return y
                return 0
            end
        end

        @test abs(1 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 10

        free_particle!.(eachrow(phasepoints), 10)
        @test abs(1 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 20

        free_particle!.(eachrow(phasepoints), 100)
        @test abs(1 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 100

        free_particle!.(eachrow(phasepoints), 1000)
        @test abs(1 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 1_000
    end

    @testset "$N Points on Quarter Circle in $(D)D " for D in [4, 12], N in [234, 5678]
        Ω(z, t, p) = CanonicalSymplecticMatrix(D)
        pinv = SecondPoincareInvariant{Float64}(Ω, D, N, Val(false))

        parampoints = getpoints(pinv)

        phasepoints = getpoints(pinv) do r, θ
            ntuple(D) do i
                i == 1 && return r * cos(θ * π/2)
                i == D ÷ 2 + 1 && return r * sin(θ * π/2)
                return 0
            end
        end

        @test abs(π/4 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 10

        free_particle!.(eachrow(phasepoints), 10)
        @test abs(π/4 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 20

        free_particle!.(eachrow(phasepoints), 100)
        @test abs(π/4 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 200

        free_particle!.(eachrow(phasepoints), 1000)
        @test abs(π/4 - compute!(pinv, phasepoints, 0, nothing)) / eps() < 2_000
    end
end
