# Known issues

## KI-1 · scope · `test/Project.toml` bounds three dependencies it already had

`test/Project.toml` gives `ChebyshevTransforms = "0.0.1"`, `GeometricEquations = "0.21"` and
`GeometricSolutions = "0.6"` the root bounds. These three were in `test/Project.toml` before the
test migration, with no bound. The plan's compat rule covers only the dependencies that the
migration adds or moves (here LinearAlgebra and StaticArrays). The bounds equal the root
`[compat]`, so no resolve changes. Found by both critics of the migration.

## KI-2 · not verified · the quality files on Julia 1.10

`test/quality/aqua.jl` and `test/quality/doctests.jl` were not run on Julia 1.10 locally: the
sandbox blocked an artifact download. The CI `min` jobs are the check. `quality/doctests.jl` runs
on every matrix entry; the doctested outputs are integer matrix displays, so a version-dependent
output is unlikely.
