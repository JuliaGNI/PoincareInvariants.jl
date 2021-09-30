@safetestset "ChebyshevImplementation" begin include("test_ChebyshevImplementation.jl") end

# @safetestset "SecondPoincareInvariant Interface" begin include("test_interface.jl") end

# @safetestset "Free Particle" begin
#     using PoincareInvariants

#     function free_particle!(state, δt)
#         mid = length(state) ÷ 2
#         state[1:mid] += state[mid+1:end] .* δt
#     end

#     @testset "Free Particle in $(D)D on 0..1 × 0..1" for D in [2, 6]
#         Ω = CanonicalSymplecticMatrix(D)
#         pinv = SecondPoincareInvariant{D, Float64}(Ω, 500)

#         parampoints = getpoints(pinv)

#         phasepoints = ones(Float64, pinv.N, D)
#         phasepoints[:, 1] .= parampoints[:, 1]
#         phasepoints[:, D ÷ 2 + 1] .= parampoints[:, 2]

#         @test compute(pinv, phasepoints) ≈ 1 atol=1e-10

#         free_particle!.(eachrow(phasepoints), 10)
#         @test compute(pinv, phasepoints) ≈ 1 atol=1e-10

#         free_particle!.(eachrow(phasepoints), 100)
#         @test compute(pinv, phasepoints) ≈ 1 atol=1e-10

#         free_particle!.(eachrow(phasepoints), 1000)
#         @test compute(pinv, phasepoints) ≈ 1 atol=1e-10
#     end

#     @testset "Free Particle in $(D)D on Quarter Circle" for D in [4, 12]
#         Ω = CanonicalSymplecticMatrix(D)
#         pinv = SecondPoincareInvariant{D, Float64}(Ω, 500)

#         parampoints = getpoints(pinv)

#         phasepoints = ones(Float64, pinv.N, D)
#         phasepoints[:, 2] .= map(eachrow(parampoints)) do row
#             row[1] * cos(row[2] * π/2)
#         end
#         phasepoints[:, D ÷ 2 + 2] .= map(eachrow(parampoints)) do row
#             row[1] * sin(row[2] * π/2)
#         end

#         @test compute(pinv, phasepoints) ≈ π/4 atol=1e-10

#         free_particle!.(eachrow(phasepoints), 10)
#         @test compute(pinv, phasepoints) ≈ π/4 atol=1e-10

#         free_particle!.(eachrow(phasepoints), 100)
#         @test compute(pinv, phasepoints) ≈ π/4 atol=1e-10

#         free_particle!.(eachrow(phasepoints), 1000)
#         @test compute(pinv, phasepoints) ≈ π/4 atol=1e-10
#     end
# end
