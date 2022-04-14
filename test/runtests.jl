using SafeTestsets, Test

@safetestset "CanonicalSymplecticForms" begin include("test_CanonicalSymplecticForms.jl") end

@safetestset "Plan Unit Tests" begin
    @safetestset "FirstFinDiffPlans" begin include("test_FirstFinDiffPlans.jl") end
    @safetestset "FirstFourierPlans" begin include("test_FirstFourierPlans.jl") end

    @safetestset "SecondChebyshevPlans" begin include("test_SecondChebyshevPlans.jl") end
    @safetestset "SecondFinDiffPlans" begin include("test_SecondFinDiffPlans.jl") end
end

@safetestset "PoincareInvariants Function Tests" begin
    include("test_PoincareInvariants.jl")
end
