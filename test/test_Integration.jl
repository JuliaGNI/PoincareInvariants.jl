using PoincareInvariants
using GeometricEquations
using GeometricIntegrators
using GeometricProblems
using StaticArrays
using Test

const HarmonicOscillator     = GeometricProblems.HarmonicOscillator
const Pendulum               = GeometricProblems.Pendulum
const MasslessChargedParticleSingular = GeometricProblems.MasslessChargedParticleSingular
const LotkaVolterra2dSingular = GeometricProblems.LotkaVolterra2dSingular

# maximum deviation of a time series of invariant values from its initial value
maxdev(Is) = maximum(abs, Is .- Is[1])

## Canonical systems ###########################################################
# HarmonicOscillator and Pendulum are canonical Hamiltonian (q, p) systems. We
# integrate their HODE form with the symplectic PartitionedGauss method, which
# preserves the canonical first (∮ p dq) and second (∫ dq ∧ dp) invariants.

@testset "HarmonicOscillator (canonical)" begin
    prob = HarmonicOscillator.hodeproblem([0.0], [0.0]; timespan = (0.0, 1.0), timestep = 0.1)

    # first invariant: a circle of radius r around the origin encloses area π r²
    r = 0.5
    pi1 = CanonicalFirstPI{Float64, 2}(500)
    ens1 = PIEnsembleProblem(prob, pi1, ϕ -> (r * sinpi(2ϕ), r * cospi(2ϕ)))
    sol1 = integrate(ens1, PartitionedGauss(1))
    Is1 = compute!(pi1, sol1)
    @test Is1[1] ≈ π * r^2 atol = 1e-12
    @test maxdev(Is1) < 1e-11

    # second invariant: a unit square centred at the origin has area one
    pi2 = CanonicalSecondPI{Float64, 2}(2_000)
    ens2 = PIEnsembleProblem(prob, pi2, (x, y) -> (x - 0.5, y - 0.5))
    sol2 = integrate(ens2, PartitionedGauss(1))
    Is2 = compute!(pi2, sol2)
    @test Is2[1] ≈ 1.0 atol = 1e-12
    @test maxdev(Is2) < 1e-11
end

@testset "Pendulum (canonical)" begin
    prob = Pendulum.hodeproblem([0.0], [0.0]; timespan = (0.0, 1.0), timestep = 0.1)

    r = 0.5
    pi1 = CanonicalFirstPI{Float64, 2}(500)
    ens1 = PIEnsembleProblem(prob, pi1, ϕ -> (r * sinpi(2ϕ), r * cospi(2ϕ)))
    sol1 = integrate(ens1, PartitionedGauss(1))
    Is1 = compute!(pi1, sol1)
    @test Is1[1] ≈ π * r^2 atol = 1e-12
    @test maxdev(Is1) < 1e-11

    pi2 = CanonicalSecondPI{Float64, 2}(2_000)
    ens2 = PIEnsembleProblem(prob, pi2, (x, y) -> (x - 0.5, y - 0.5))
    sol2 = integrate(ens2, PartitionedGauss(1))
    Is2 = compute!(pi2, sol2)
    @test Is2[1] ≈ 1.0 atol = 1e-12
    @test maxdev(Is2) < 1e-11
end

## Noncanonical systems ########################################################
# MasslessChargedParticleSingular and LotkaVolterra2dSingular carry a state-dependent
# (noncanonical) symplectic structure. We integrate their degenerate IODE/LODE form
# with a Degenerate Variational Runge-Kutta method (DVRK), which is designed for
# noncanonical symplectic systems and preserves the noncanonical one-form ϑ and
# two-form ω = dϑ (unlike the DVI methods, which preserve a discrete canonical
# structure). Both use a "singular" gauge whose one-form is exactly the symplectic
# potential DVRK expects, so the invariants are preserved to machine accuracy. The
# invariant forms are reused directly from the GeometricProblems modules.

@testset "MasslessChargedParticleSingular (noncanonical)" begin
    prob = MasslessChargedParticleSingular.lodeproblem([1.0, 1.0]; timespan = (0.0, 1.0), timestep = 0.05)
    par  = parameters(prob)

    # the in-place one-form ϑ and two-form ω = dϑ of the problem match the invariant form
    # convention `form(out, t, z, p)` exactly, so we pass them straight to FirstPI / SecondPI
    # with no wrapper — identical setup to LotkaVolterra2dSingular below
    fs = functions(prob)

    q₀ = SVector(1.0, 1.0)
    ρ = 0.2

    pi1 = FirstPI{Float64, 2}(fs.ϑ, 500)
    ens1 = PIEnsembleProblem(prob, pi1, ϕ -> q₀ .+ ρ .* (cospi(2ϕ), sinpi(2ϕ)))
    sol1 = integrate(ens1, DVRK(Gauss(2)))
    Is1 = compute!(pi1, sol1, par)
    @test maxdev(Is1) < 1e-13

    pi2 = SecondPI{Float64, 2}(fs.ω, 2_000)
    ens2 = PIEnsembleProblem(prob, pi2, (x, y) -> q₀ .+ ρ .* (2x - 1, 2y - 1))
    sol2 = integrate(ens2, DVRK(Gauss(2)))
    Is2 = compute!(pi2, sol2, par)
    @test maxdev(Is2) < 1e-13
end

@testset "LotkaVolterra2dSingular (noncanonical)" begin
    # the "singular" Lagrangian carries the form of the symplectic potential that DVRK
    # expects, so the noncanonical invariants are preserved to machine accuracy
    prob = LotkaVolterra2dSingular.lodeproblem([2.0, 1.0]; timespan = (0.0, 1.0), timestep = 0.05)
    par  = parameters(prob)

    # the problem's in-place one-form ϑ and two-form ω = dϑ follow the `form(out, t, z, p)`
    # convention, so they are passed directly as the invariant forms
    fs = functions(prob)

    q₀ = SVector(2.0, 1.0)
    ρ = 0.2

    pi1 = FirstPI{Float64, 2}(fs.ϑ, 500)
    ens1 = PIEnsembleProblem(prob, pi1, ϕ -> q₀ .+ ρ .* (cospi(2ϕ), sinpi(2ϕ)))
    sol1 = integrate(ens1, DVRK(Gauss(2)))
    Is1 = compute!(pi1, sol1, par)
    @test maxdev(Is1) < 1e-11

    pi2 = SecondPI{Float64, 2}(fs.ω, 2_000)
    ens2 = PIEnsembleProblem(prob, pi2, (x, y) -> q₀ .+ ρ .* (2x - 1, 2y - 1))
    sol2 = integrate(ens2, DVRK(Gauss(2)))
    Is2 = compute!(pi2, sol2, par)
    @test maxdev(Is2) < 1e-11
end
