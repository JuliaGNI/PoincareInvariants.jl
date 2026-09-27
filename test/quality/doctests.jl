using Documenter
using PoincareInvariants

# Documenter evaluates the `@meta` block of a manual page in `Main`, and this file runs in a
# module of its own.
@eval Main import PoincareInvariants

# The same setup as `docs/make.jl` and the CI Doctests job.
include(joinpath(pkgdir(PoincareInvariants), "docs", "doctestsetup.jl"))

doctest(PoincareInvariants)
