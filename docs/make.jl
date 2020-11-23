using Documenter, PoincareInvariants

makedocs(
    sitename = "PoincareInvariants.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = ["Home" => "index.md",
             "Poincaré Integral Invariants" => "poincare_invariants.md",
             "1st Poincaré Invariant" => "poincare_invariant_1st.md",
             "2nd Poincaré Invariant" => "poincare_invariant_2nd.md",
             ]
)

deploydocs(
    repo   = "github.com/JuliaGNI/PoincareInvariants.jl"
)
