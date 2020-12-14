
include("poincare_invariant_2nd_module.jl")

using BenchmarkTools
using GeometricIntegrators
using Test

using .PoincareInvariant2ndTest

set_config(:verbosity, 1)

const nx = 141
const ny = 142


println("\nSecond Canonical Euler-Poincaré Integral Invariant (ApproxFun):")
@btime compute_canonical_invariant_approxfun(nx, ny)

println("\nSecond Euler-Poincaré Integral Invariant (ApproxFun):")
@btime compute_invariant_approxfun(nx, ny)

println("\nSecond Euler-Poincaré Integral Invariant (OrthogonalPolynomials):")
@btime compute_invariant_opq(nx, ny)

println("\nSecond Euler-Poincaré Integral Invariant (Trapezoidal):")
@btime compute_invariant_trapezoidal(nx, ny)
