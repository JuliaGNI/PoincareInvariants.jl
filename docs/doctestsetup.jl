# What the doctests of this package need in scope.
#
# Included both by `docs/make.jl` and by the `doctest` job of `.github/workflows/CI.yml`, so that a
# documentation build and a doctest run cannot disagree. One definition, two callers: the CI
# workflow is byte-identical in every repository and cannot carry per-package knowledge, and a
# second copy of this list is exactly the thing that goes stale.
#
# Two calls rather than one recursive one: the doctests of `CanonicalSymplecticForms` are written
# against that submodule's own exports, not against the parent's, so `recursive = true` from
# `PoincareInvariants` would put the wrong names in scope for them.

using Documenter: DocMeta

using PoincareInvariants

DocMeta.setdocmeta!(
    PoincareInvariants,
    :DocTestSetup,
    quote
        using PoincareInvariants
    end
)

DocMeta.setdocmeta!(
    PoincareInvariants.CanonicalSymplecticForms,
    :DocTestSetup,
    quote
        using PoincareInvariants.CanonicalSymplecticForms
    end
)
