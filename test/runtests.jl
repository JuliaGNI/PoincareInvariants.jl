using SafeTestsets, Test

@safetestset "CanoincalSymplecticMatrices" begin include("test_CanonicalSymplecticMatrices.jl") end
@safetestset "FirstPoincareInvariants" begin include("test_FirstPoincareInvariants.jl") end
@safetestset "SecondPoincareInvariants" begin include("test_SecondPoincareInvariants.jl") end

# include("padua.jl")
# include("second_poincare_invariant.jl")
