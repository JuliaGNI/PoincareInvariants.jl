using Documenter, PoincareInvariants

makedocs(
    sitename = "PoincareInvariants.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home" => "index.md",
        "First Poincaré Invariants" => "first_poincare_invariants.md",
        "Second Poincaré Invariants" => "second_poincare_invariants.md",
        "Canonical Symplectic Structures" => "canonical_symplectic_structures.md",
        "Reference" => "reference.md"
    ]
)

deploydocs(
    repo   = "github.com/JuliaGNI/PoincareInvariants.jl"
)
