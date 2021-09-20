using SafeTestsets, Test

@safetestset "FirstPoincareInvariants" begin
    include("test_FirstPoincareInvariants.jl")
end

@safetestset "SecondPoincareInvariants" begin
    include("SecondPoincareInvariants/test_SecondPoincareInvariants.jl")
end

# include("padua.jl")
# include("second_poincare_invariant.jl")
