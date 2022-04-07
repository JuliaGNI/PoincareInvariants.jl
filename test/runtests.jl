using SafeTestsets, Test

@safetestset "CanonicalSymplecticStructures" begin
    include("test_CanonicalSymplecticStructures.jl")
end

@safetestset "Plan Unit Tests" begin
    @safetestset "FirstFinDiffPlans" begin include("test_FirstFinDiffPlans.jl") end

    @safetestset "SecondChebyshevPlans" begin include("test_SecondChebyshevPlans.jl") end
    @safetestset "SecondFinDiffPlans" begin include("test_SecondFinDiffPlans.jl") end
end

@safetestset "PoincareInvariants Function Tests" begin
    include("test_PoincareInvariants.jl")
end
