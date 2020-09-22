
using GeometricIntegrators
using SafeTestsets

set_config(:verbosity, 2)

@safetestset "1st Poincaré Invariant Unit Tests                                               " begin include("poincare_invariant_1st_unittests.jl") end
@safetestset "1st Poincaré Invariant Function Tests                                           " begin include("poincare_invariant_1st_tests.jl") end
@safetestset "2nd Poincaré Invariant Unit Tests                                               " begin include("poincare_invariant_2nd_unittests.jl") end
@safetestset "2nd Poincaré Invariant Function Tests                                           " begin include("poincare_invariant_2nd_tests.jl") end
