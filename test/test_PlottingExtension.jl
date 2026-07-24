using PoincareInvariants
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using CairoMakie   # loads Makie, activating PoincareInvariantsMakieExt
using Test

# small canonical problem for a quick smoke test of the plotting extension
probh = hodeproblem([0.0], [0.0]; timespan = (0.0, 1.0), timestep = 0.2)

init1 = ϕ -> (0.5 * sinpi(2ϕ), 0.5 * cospi(2ϕ))
init2 = (x, y) -> (x - 0.5, y - 0.5)

pi1 = CanonicalFirstPI{Float64, 2}(50)
sol1 = integrate(PIEnsembleProblem(probh, pi1, init1), PartitionedGauss(1))

pi2 = CanonicalSecondPI{Float64, 2}(100)
sol2 = integrate(PIEnsembleProblem(probh, pi2, init2), PartitionedGauss(1))

grid = CanonicalSecondPI{Float64, 2}((7, 7), SecondFinDiffPlan)
solg = integrate(PIEnsembleProblem(probh, grid, init2), PartitionedGauss(1))

@testset "invariant error plots" begin
    @test plot_invariant(pi1, sol1) isa Figure
    @test plot_invariant(pi2, sol2; title = "second") isa Figure
    @test plot_invariant(pi1, "a" => sol1, "b" => sol1) isa Figure
end

@testset "advected loop and surface" begin
    @test plot_loop(sol1) isa Figure
    @test plot_surface(grid, solg) isa Figure
end
