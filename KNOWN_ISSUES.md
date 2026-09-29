# Known issues

### K2 · `test/quality/aqua.jl` and `test/quality/doctests.jl` have no local run on Julia 1.10.

- **location:** `test/quality/aqua.jl`
- **evidence:** the sandbox blocked an artifact download. The CI `min` jobs are the check.
  `quality/doctests.jl` runs on every matrix entry; the doctested outputs are integer matrix
  displays (`src/CanonicalSymplecticForms.jl:41`), so a version-dependent output is unlikely.
- **kind:** not verified
- **found:** 2026-09-26
