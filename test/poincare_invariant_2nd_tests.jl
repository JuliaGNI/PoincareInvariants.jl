
include("poincare_invariant_2nd_module.jl")

using Test

using .PoincareInvariant2ndTest

const nx = 141
const ny = 142


### compute analytical invariant ###

Ia = compute_invariant_analytical()


### compute and check canonical invariant ###

Jn = compute_canonical_invariant_approxfun(nx, ny)

@test Jn ≈ Ia atol=2eps()


### compute and check noncanonical invariant ###

In, Jn = compute_invariant_approxfun(nx, ny)

@test In ≈ Ia atol=2eps()
@test Jn ≈ Ia atol=2eps()


In, Jn = compute_invariant_opq(nx, ny)

@test In ≈ Ia atol=1E-1
#@test Jn ≈ Ia atol=1E-1


In, Jn = compute_invariant_trapezoidal(nx, ny)

@test In ≈ Ia atol=1E-2
@test Jn ≈ Ia atol=1E-2
