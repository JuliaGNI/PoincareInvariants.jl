# Changelog

All notable changes to `PoincareInvariants.jl` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project follows [Semantic Versioning](https://semver.org/) — with the usual 0.x convention that
a **minor** version bump may contain breaking changes.

## [0.5.0] - 2026-07-24

This release moves `PoincareInvariants` onto the [JuliaGNI](https://github.com/JuliaGNI)
ecosystem and reworks the differential-form interface. It contains several breaking changes.

### Breaking

- **Migrated from the SciML ecosystem to JuliaGNI.** Dynamics are now described with
  [GeometricEquations.jl](https://github.com/JuliaGNI/GeometricEquations.jl) problems
  (`ODEProblem`, `HODEProblem`, `PODEProblem`, `IODEProblem`, `LODEProblem`) and integrated with
  [GeometricIntegrators.jl](https://github.com/JuliaGNI/GeometricIntegrators.jl); solutions are
  `GeometricSolutions.EnsembleSolution`s. The `SciMLBase` and `CommonSolve` dependencies were
  removed and replaced by `GeometricEquations` and `GeometricSolutions`.
- **Integration entry point changed.** The old `solve(::PIEnsembleProblem, alg, ensemblealg)` is
  removed. Build the ensemble and integrate with GeometricIntegrators instead:
  `sol = integrate(PIEnsembleProblem(prob, pinv, init), method)`, then `compute!(pinv, sol)`.
- **`PIEnsembleProblem` argument order changed** from `PIEnsembleProblem(init, prob, pinv)` to
  `PIEnsembleProblem(prob, pinv, init)`.
- **Differential forms are now in-place.** A form is supplied as `form(out, t, z, p)` — it writes
  its value at the phase-space point `z` into the preallocated output `out` (a length-`D` vector
  for a one-form `ϑ`, a `D×D` matrix for a two-form `ω`) — rather than the previous
  value-returning `form(z, t, p)`. The argument order matches the forms exported by
  GeometricEquations/GeometricProblems (`ϑ(Θ, t, q, params)`, `ω(Ω, t, q, params)`), so a
  problem's own forms can be passed directly with no wrapper. A constant `ω::AbstractMatrix`
  (e.g. a `CanonicalSymplecticMatrix`) is still accepted and used as-is.
- **Canonical-form API changed.** Removed the value-returning `canonical_one_form` and
  `canonical_two_form`, and the lazy `CanonicalSymplecticVector`; added the in-place
  `canonical_one_form!`. `CanonicalSymplecticMatrix` is unchanged and remains the constant
  canonical two-form.
- **Julia 1.10 or newer is now required** (previously 1.6).

### Added

- Plotting via a Makie package extension: `plot_invariant`, `plot_loop` and `plot_surface`
  become available once a Makie backend (e.g. `CairoMakie`) is loaded.
- Documentation on the theory and implementation, a rewritten tutorial covering the in-place
  form interface with worked canonical and noncanonical examples, and a convergence study.
- Noncanonical examples built on `GeometricProblems` singular problems
  (`MasslessChargedParticleSingular`, `LotkaVolterra2dSingular`), whose one- and two-forms are
  taken directly from the `LODEProblem` via `functions(prob)`.

### Performance

- Forms are evaluated in place into a single reused `StaticArrays` buffer (`MVector`/`MMatrix`),
  so computing an invariant no longer allocates per sample point.

[0.5.0]: https://github.com/JuliaGNI/PoincareInvariants.jl/compare/v0.4.0...v0.5.0
