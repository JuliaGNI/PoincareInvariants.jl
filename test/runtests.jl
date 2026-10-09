using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "CanonicalSymplecticForms" include("CanonicalSymplecticForms.jl")
    @safetestset "FirstFinDiffPlans" include("FirstFinDiffPlans.jl")
    @safetestset "FirstFourierPlans" include("FirstFourierPlans.jl")
    @safetestset "SecondChebyshevPlans" include("SecondChebyshevPlans.jl")
    @safetestset "SecondFinDiffPlans" include("SecondFinDiffPlans.jl")
    @safetestset "PoincareInvariants" include("PoincareInvariants.jl")
    @safetestset "Integration with GeometricIntegrators" include("integration/geometric_integrators.jl")
    @safetestset "Plotting extension" include("integration/makie_extension.jl")
end
if "doctests" in GROUPS
    @safetestset "Doctests" include("quality/doctests.jl")
end
