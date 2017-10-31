using Documenter, PoincareInvariants

makedocs(
    sitename = "PoincareInvariants.jl",
    format = :html,
    pages = ["Home" => "index.md",
             "Poincaré Invariants"    => "poincare_invariants.md",
             "1st Poincaré Invariant" => "poincare_invariant_1st.md",
             "2nd Poincaré Invariant" => "poincare_invariant_2nd.md",
             ]
)

deploydocs(
    repo   = "github.com/DDMGNI/PoincareInvariants.jl.git",
    target = "build",
    julia  = "0.6",
    osname = "linux",
    deps   = nothing,
    make   = nothing)

# note: julia version must be the same as in .travis.yml, that is both must be
#       set to either 0.6 or release, but not one to 0.6 and one to release.
