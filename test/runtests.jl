using SafeTestsets, Test

@safetestset "CanonicalSymplecticStructures" begin
    include("test_CanonicalSymplecticStructures.jl")
end

@safetestset "FirstPoincareInvariants" begin
    include("FirstPoincareInvariants/test_FirstPoincareInvariants.jl")
end

@safetestset "SecondPoincareInvariants" begin
    include("SecondPoincareInvariants/test_SecondPoincareInvariants.jl")
end
