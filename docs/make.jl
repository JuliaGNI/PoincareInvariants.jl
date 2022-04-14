using Documenter, PoincareInvariants

DocMeta.setdocmeta!(
    PoincareInvariants,
    :DocTestSetup,
    quote
        using PoincareInvariants
    end
)

DocMeta.setdocmeta!(
    PoincareInvariants.SecondChebyshevPlans.PaduaTransforms,
    :DocTestSetup,
    quote
        using PoincareInvariants.SecondChebyshevPlans.PaduaTransforms
    end
)

DocMeta.setdocmeta!(
    PoincareInvariants.CanonicalSymplecticForms,
    :DocTestSetup,
    quote
        using PoincareInvariants.CanonicalSymplecticForms
    end
)

makedocs(
    sitename = "PoincareInvariants.jl",
    modules=[PoincareInvariants],
    pages = [
        "Home" => "index.md",
        "First Poincaré Invariants" => "first_poincare_invariants.md",
        "Second Poincaré Invariants" => [
            "Guide" => "second_poincare_invariants.md",
            "Chebyshev Implementation" => "chebyshev_implementation.md",
            "Padua Transforms" => "padua_transforms.md"
        ],
        "Canonical Symplectic Structures" => "canonical_symplectic_structures.md",
        "Reference" => "reference.md"
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo   = "github.com/JuliaGNI/PoincareInvariants.jl"
)
