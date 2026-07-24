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
        "Tutorial" => "tutorial.md",
        "Theory" => "theory.md",
        "Implementation" => "implementation.md",
        "Convergence" => "convergence.md",
        "Examples" => [
            "Harmonic Oscillator" => "examples/harmonic_oscillator.md",
            "Pendulum" => "examples/pendulum.md",
            "Massless Charged Particle" => "examples/massless_charged_particle.md",
            "Lotka-Volterra" => "examples/lotka_volterra.md",
        ],
        "Reference" => "reference.md"
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(repo = "github.com/JuliaGNI/PoincareInvariants.jl")
