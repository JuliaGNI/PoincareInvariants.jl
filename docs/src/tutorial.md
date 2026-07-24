# Computing Invariants with GeometricIntegrators.jl

`PoincareInvariants.jl` is built on the [JuliaGNI](https://github.com/JuliaGNI) ecosystem.
Problems are defined with [GeometricEquations.jl](https://github.com/JuliaGNI/GeometricEquations.jl),
integrated with [GeometricIntegrators.jl](https://github.com/JuliaGNI/GeometricIntegrators.jl),
and a large collection of ready-made example systems is provided by
[GeometricProblems.jl](https://github.com/JuliaGNI/GeometricProblems.jl).

The general workflow for computing a Poincaré integral invariant of a dynamical system is:

1. **Choose a problem.** Any geometric `EquationProblem` — an `ODEProblem`, `HODEProblem`,
   `PODEProblem`, `IODEProblem` or `LODEProblem` — specifies the dynamics, time span and time
   step.
2. **Set up an invariant.** Create a `FirstPI`/`SecondPI` (or the canonical variants
   `CanonicalFirstPI`/`CanonicalSecondPI`) for the required phase space dimension and
   differential form.
3. **Build an ensemble.** [`PIEnsembleProblem`](@ref) samples the initial curve or surface
   with [`getpoints`](@ref) and turns each sampled phase space point into the initial
   condition of one member of a `GeometricEquations.EnsembleProblem`.
4. **Integrate.** Call `integrate(ensemble, method)` from GeometricIntegrators, which returns
   an `EnsembleSolution`.
5. **Compute.** Pass the solution to [`compute!`](@ref) to obtain the invariant at every saved
   time step.

## Building the ensemble

[`PIEnsembleProblem`](@ref)`(prob, pinv, init)` takes the base problem `prob`, the invariant setup
object `pinv`, and the curve or surface parameterisation `init`. Each sampled point becomes
one initial condition:

- for an `ODEProblem` the point is the position `q`;
- for a `HODEProblem`/`PODEProblem` the point is split into its position and momentum halves
  `(q, p)`;
- for an `IODEProblem`/`LODEProblem` the position is taken from the point and the momentum is
  initialised from the equation's one-form as `p₀ = ϑ(t₀, q₀)`.

The parameterisation `init` is itself a function that **returns** the sampled point (this is
separate from the differential form discussed below):

- for the **first** invariant the loop is parameterised by one variable, `ϕ -> (…)`, and must
  return `D` values (one per phase space coordinate);
- for the **second** invariant the surface is parameterised by two variables, `(x, y) -> (…)`,
  again returning `D` values.

## Differential forms

The invariant is the integral of a differential form over the advected curve or surface. How
you supply that form depends on whether the system is canonical or noncanonical.

- **Canonical** Hamiltonian systems (`HODEProblem`/`PODEProblem` in coordinates `(q, p)`) use
  the *canonical* symplectic form. Use `CanonicalFirstPI` / `CanonicalSecondPI`,
  which supply that form internally — you pass **no** form argument.
- **Noncanonical** symplectic systems carry a state-dependent one-form `ϑ` (and two-form
  `ω = dϑ`). Use `FirstPI` / `SecondPI` and pass the form explicitly.

A user-supplied form is an **in-place** function following the convention

```julia
form(out, t, z, p)
```

where `out` is the preallocated output the form fills — a length-`D` vector for a one-form
`ϑ`, a `D×D` matrix for a two-form `ω` — `z` is one phase space point, `t` the time and `p`
the parameters. This argument order is exactly the one `GeometricEquations`/`GeometricProblems`
use for their forms (`ϑ(Θ, t, q, params)`, `ω(Ω, t, q, params)`), so **a problem's own forms
can be passed directly**, with no wrapper (see the noncanonical example below). Internally the
form is only ever called with a preallocated buffer that is reused across sample points, so
evaluating it allocates nothing.

## Integrating and computing

Because the number of ensemble members equals the number of sampled points, no trajectory
count needs to be supplied — simply integrate and compute. Here is the full workflow for the
first invariant of a **canonical** system, the harmonic oscillator:

```@example tutorial_canonical
using PoincareInvariants
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator

prob = hodeproblem([0.0], [0.0]; timespan = (0.0, 1.0), timestep = 0.1)

# a circle of radius r around the origin encloses the area π r²
r = 0.5
pinv = CanonicalFirstPI{Float64, 2}(500)
ens  = PIEnsembleProblem(prob, pinv, ϕ -> (r * sinpi(2ϕ), r * cospi(2ϕ)))

sol = integrate(ens, PartitionedGauss(1))
Is  = compute!(pinv, sol)

Is[1], maximum(abs, Is .- Is[1])   # initial value ≈ π r², and its drift
```

[`compute!`](@ref) returns a `Vector` holding one invariant value per saved time step. A
symplectic integrator keeps that value essentially constant, so its drift from `Is[1]` (here
`maximum(abs, Is .- Is[1])`) is the diagnostic of interest — not the individual values.

The second invariant works the same way with `CanonicalSecondPI` and a surface
parameterisation `(x, y) -> (…)`; see the [Harmonic Oscillator](examples/harmonic_oscillator.md)
and [Pendulum](examples/pendulum.md) pages.

## A noncanonical example

For a **noncanonical** system we integrate the degenerate `LODEProblem` (or `IODEProblem`) and
supply the problem's own one-form/two-form. Both are obtained from the problem with
`functions(prob)` and passed directly — `fs.ϑ` and `fs.ω` already have the in-place
`form(out, t, z, p)` signature:

```@example tutorial_noncanonical
using PoincareInvariants
using GeometricIntegrators
using GeometricProblems.LotkaVolterra2dSingular

prob = lodeproblem([2.0, 1.0]; timespan = (0.0, 1.0), timestep = 0.02)
par  = parameters(prob)   # physical parameters carried by the problem
fs   = functions(prob)    # the problem's own ϑ, ω, …

q₀ = (2.0, 1.0)
ρ  = 0.2

# first invariant: advect a small loop around q₀ and pass fs.ϑ directly
pi1  = FirstPI{Float64, 2}(fs.ϑ, 500)
ens1 = PIEnsembleProblem(prob, pi1, ϕ -> q₀ .+ ρ .* (cospi(2ϕ), sinpi(2ϕ)))
sol1 = integrate(ens1, DVRK(Gauss(2)))
Is1  = compute!(pi1, sol1, par)

# second invariant: advect a small square and pass fs.ω directly
pi2  = SecondPI{Float64, 2}(fs.ω, 2_000)
ens2 = PIEnsembleProblem(prob, pi2, (x, y) -> q₀ .+ ρ .* (2x - 1, 2y - 1))
sol2 = integrate(ens2, DVRK(Gauss(2)))
Is2  = compute!(pi2, sol2, par)

maximum(abs, Is1 .- Is1[1]), maximum(abs, Is2 .- Is2[1])   # both ≈ machine precision
```

Note that `par` is passed to [`compute!`](@ref) because the forms depend on the physical
parameters; omitting it for a parameter-dependent form is a common mistake (see below). The
fuller [Lotka-Volterra](examples/lotka_volterra.md) and
[Massless Charged Particle](examples/massless_charged_particle.md) pages carry this through to
plots of the invariant over time.

## Choosing an invariant and a plan

Each invariant can be computed with more than one numerical plan (see
[Implementation](implementation.md) for the details):

| invariant | default plan | alternative |
|-----------|--------------|-------------|
| first (`FirstPI`)  | `FirstFourierPlan` (spectral, periodic loop) | `FirstFinDiffPlan` (finite differences) |
| second (`SecondPI`)| `SecondChebyshevPlan` (spectral, Padua points) | `SecondFinDiffPlan` (finite differences) |

The plan is the optional last constructor argument, e.g. `FirstPI{Float64, 2}(ϑ, 500, FirstFinDiffPlan)`.
The default spectral plans converge exponentially for smooth data and are almost always the
right choice. `SecondFinDiffPlan` additionally accepts a non-square grid as a tuple `(Nx, Ny)`
instead of a single point count `N`.

## Choosing an integrator

The invariant is only conserved if the integrator preserves the relevant symplectic
structure.

- For **canonical** Hamiltonian systems (`HODEProblem`/`PODEProblem`) a symplectic partitioned
  method such as `PartitionedGauss` conserves the canonical invariants to machine precision.
- For **noncanonical** symplectic systems the degenerate `IODEProblem`/`LODEProblem`
  formulation should be integrated with a **Degenerate Variational Runge-Kutta** method,
  e.g. `DVRK(Gauss(2))`, which preserves the noncanonical symplectic structure. (The plain
  DVI methods preserve a discrete *canonical* structure and are not suitable here.)

## Common pitfalls

- **Forgetting the parameters.** If the differential form depends on physical parameters, call
  `compute!(pinv, sol, p)` (and `plot_invariant(pinv, sol; p = …)`). Omitting `p` leaves the
  form without its parameters and gives wrong values. Canonical invariants take no parameters,
  so `compute!(pinv, sol)` is enough.
- **Mismatched dimension.** The `D` in `FirstPI{T, D}` / `CanonicalFirstPI{T, D}` is the phase
  space dimension and must match the problem; the parameterisation must return `D` values per
  point.
- **A non-symplectic (or wrong) integrator.** With a structure-agnostic method the invariant
  drifts even though the code runs without error. Match the integrator to the structure as
  above; judge conservation by the drift of `compute!`'s output from its initial value.
- **Reading individual values instead of the drift.** `compute!` returns the invariant at every
  saved time step. What matters is how little it changes over time, i.e. `Is .- Is[1]`.

The following guides work through a canonical example
([Harmonic Oscillator](examples/harmonic_oscillator.md), [Pendulum](examples/pendulum.md)) and a
noncanonical example ([Massless Charged Particle](examples/massless_charged_particle.md),
[Lotka-Volterra](examples/lotka_volterra.md)) in detail.
