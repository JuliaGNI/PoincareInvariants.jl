
include("poincare_invariant_2nd_module.jl")

using Requires
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


### compute and check noncanonical invariant computed with Orthogonal Polynomials ###

# @require ClassicalOrthogonalPolynomials = "b30e2e7b-c4ee-47da-9d5f-2c5c27239acd" begin
#     @require FastTransforms = "057dd010-8810-581a-b7be-e3fc3b93f78c" begin
        In, Jn = compute_invariant_cop(nx, ny)

        @test In ≈ Ia atol=1E-2
        # @test Jn ≈ Ia atol=1E-2
#     end
# end


### compute and check noncanonical invariant computed with Finite Differences ###

In, Jn = compute_invariant_trapezoidal(nx, ny)

@test In ≈ Ia atol=1E-2
@test Jn ≈ Ia atol=1E-2
