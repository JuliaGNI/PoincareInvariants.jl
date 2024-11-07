using Documenter, PoincareInvariants

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

makedocs(
    sitename = "PoincareInvariants.jl",
    modules=[PoincareInvariants],
    warnonly = Documenter.except(:autodocs_block, :cross_references, :docs_block, :doctest, :eval_block, :example_block, :footnote, :linkcheck_remotes, :linkcheck, :meta_block, :parse_error, :setup_block),
    pages = [
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Guides" => [
            "Pendulum" => "guides/pendulum.md",
            "Charged Particle" => "guides/charged_particle.md",
            "DifferentialEquations.jl" => "guides/diffeq.md",
            "Plans" => "guides/plans.md",
        ],
        "Reference" => "reference.md"
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(repo = "github.com/JuliaGNI/PoincareInvariants.jl")
