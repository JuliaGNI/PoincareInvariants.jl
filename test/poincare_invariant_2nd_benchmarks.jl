
include("poincare_invariant_2nd_module.jl")

using BenchmarkTools
using GeometricIntegrators
using Requires
using Test

using .PoincareInvariant2ndTest

set_config(:verbosity, 1)

const nx = 141
const ny = 142


println("\nSecond Canonical Euler-Poincaré Integral Invariant (ApproxFun):")
@btime compute_canonical_invariant_approxfun(nx, ny)

println("\nSecond Euler-Poincaré Integral Invariant (ApproxFun):")
@btime compute_invariant_approxfun(nx, ny)

# @require ClassicalOrthogonalPolynomials = "b30e2e7b-c4ee-47da-9d5f-2c5c27239acd" begin
#     @require FastTransforms = "057dd010-8810-581a-b7be-e3fc3b93f78c" begin
        println("\nSecond Euler-Poincaré Integral Invariant (OrthogonalPolynomials):")
        @btime compute_invariant_cop(nx, ny)
#     end
# end

println("\nSecond Euler-Poincaré Integral Invariant (Trapezoidal):")
@btime compute_invariant_trapezoidal(nx, ny)
