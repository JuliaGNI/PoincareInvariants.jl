# Known issues

### K1 · `test/Project.toml` bounds three dependencies that have no bound on `origin/main`.

- **location:** `test/Project.toml:20`
- **evidence:** lines 20–22 give `ChebyshevTransforms = "0.0.1"`, `GeometricEquations = "0.21"`
  and `GeometricSolutions = "0.6"` the root bounds. `git show origin/main:test/Project.toml` lists
  these three in `[deps]` and has no `[compat]` section. The compat rule of the test migration
  bounds only the dependencies that it adds to `test/Project.toml` or moves into it: here
  LinearAlgebra and StaticArrays. The bounds equal the root `[compat]`, so no resolve changes.
- **kind:** defect
- **found:** 2026-09-26, by both critics of the test migration

### K2 · `test/quality/aqua.jl` and `test/quality/doctests.jl` have no local run on Julia 1.10.

- **location:** `test/quality/aqua.jl`
- **evidence:** the sandbox blocked an artifact download. The CI `min` jobs are the check.
  `quality/doctests.jl` runs on every matrix entry; the doctested outputs are integer matrix
  displays (`src/CanonicalSymplecticForms.jl:41`), so a version-dependent output is unlikely.
- **kind:** not verified
- **found:** 2026-09-26

### K3 · `docs/doctestsetup.jl` names two callers and has three.

- **location:** `docs/doctestsetup.jl:3`
- **evidence:** lines 3–6 name `docs/make.jl` and the `doctest` job of `.github/workflows/CI.yml`
  as the callers ("One definition, two callers"). `test/quality/doctests.jl:9` includes the file
  too.
- **kind:** docs
- **found:** 2026-09-26
