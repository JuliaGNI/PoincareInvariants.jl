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

## Integrating and computing

Because the number of ensemble members equals the number of sampled points, no trajectory
count needs to be supplied — simply integrate and compute:

```julia
using PoincareInvariants
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator

prob = hodeproblem([0.0], [0.0]; timespan = (0.0, 1.0), timestep = 0.1)

pinv = CanonicalFirstPI{Float64, 2}(500)
ens  = PIEnsembleProblem(prob, pinv, ϕ -> (0.5 * sinpi(2ϕ), 0.5 * cospi(2ϕ)))

sol = integrate(ens, PartitionedGauss(1))
Is  = compute!(pinv, sol)
```

[`compute!`](@ref) returns a `Vector` holding one invariant value per saved time step. An
optional parameter `p` may be passed as `compute!(pinv, sol, p)`; it is handed to the
differential form together with the time (this is needed for forms that depend on physical
parameters, such as the noncanonical examples).

## Choosing an integrator

The invariant is only conserved if the integrator preserves the relevant symplectic
structure.

- For **canonical** Hamiltonian systems (`HODEProblem`/`PODEProblem`) a symplectic partitioned
  method such as `PartitionedGauss` conserves the canonical invariants to machine precision.
- For **noncanonical** symplectic systems the degenerate `IODEProblem`/`LODEProblem`
  formulation should be integrated with a **Degenerate Variational Runge-Kutta** method,
  e.g. `DVRK(Gauss(2))`, which preserves the noncanonical symplectic structure. (The plain
  DVI methods preserve a discrete *canonical* structure and are not suitable here.)

The following guides work through a canonical example
([Harmonic Oscillator](examples/harmonic_oscillator.md), [Pendulum](examples/pendulum.md)) and a noncanonical
example ([Massless Charged Particle](examples/massless_charged_particle.md),
[Lotka-Volterra](examples/lotka_volterra.md)) in detail.
