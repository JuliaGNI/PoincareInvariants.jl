# Lotka-Volterra

The two-dimensional Lotka-Volterra system is a noncanonical Hamiltonian (Poisson) system on
the positive quadrant $q = (q_1, q_2)$, $q_i > 0$, with Hamiltonian

```math
H(q) = a_1 \, q_1 + a_2 \, q_2 + b_1 \log q_1 + b_2 \log q_2 .
```

We use the `LotkaVolterra2dSingular` variant from
[GeometricProblems.jl](https://github.com/JuliaGNI/GeometricProblems.jl), whose "singular"
degenerate Lagrangian carries the one-form

```math
\vartheta(q) = \begin{pmatrix} \dfrac{\log q_2}{q_1} \\[1ex] 0 \end{pmatrix} .
```

This is exactly the form of the symplectic potential that the `DVRK` (Degenerate Variational
Runge-Kutta) integrator expects, so integrating the degenerate `LODEProblem` with `DVRK`
preserves the noncanonical symplectic structure $\omega = d\vartheta$ — and with it the
Poincaré invariants — to machine accuracy. The plotting functions `plot_loop`, `plot_surface`
and `plot_invariant` come from the package's Makie extension, activated by loading
`CairoMakie`.

```@example lotka
using PoincareInvariants
using GeometricIntegrators
using GeometricProblems.LotkaVolterra2dSingular
using StaticArrays
using CairoMakie

prob = lodeproblem([2.0, 1.0]; timespan = (0.0, 1.0), timestep = 0.02)
nothing # hide
```

## The invariant forms

The one-form $\vartheta$ and the two-form $\omega = d\vartheta$ are exported by the problem
module (`ϑ` and the in-place `ω`); we only adapt them to the `form(z, t, p)` interface.
Neither depends on the physical parameters, so `p` is ignored.

```@example lotka
oneform(z, t, p) = SVector{2}(ϑ(t, z))

function twoform(z, t, p)
    Ω = zeros(eltype(z), 2, 2)
    ω(Ω, t, z)
    SMatrix{2, 2}(Ω)
end
nothing # hide
```

## First invariant

```math
I_{1} = \oint_{\gamma} \vartheta_i(q) \, dq^i
```

We advect a circle of radius $\rho = 0.2$ around the initial point $(2, 1)$.

```@example lotka
q₀ = SVector(2.0, 1.0)
ρ  = 0.2

pi1  = FirstPI{Float64, 2}(oneform, 500)
sol1 = integrate(PIEnsembleProblem(prob, pi1, ϕ -> q₀ .+ ρ .* (cospi(2ϕ), sinpi(2ϕ))), DVRK(Gauss(2)))
nothing # hide
```

The loop is transported and deformed by the Lotka-Volterra flow:

```@example lotka
plot_loop(sol1; xlabel = "q₁", ylabel = "q₂")
```

```@example lotka
plot_invariant(pi1, sol1; title = "First Poincaré invariant")
```

## Second invariant

```math
I_{2} = \int_{S} \omega_{ij}(q) \, dq^i \, dq^j
```

We advect a square around $(2, 1)$, staying within the positive quadrant.

```@example lotka
pi2  = SecondPI{Float64, 2}(twoform, 2_000)
sol2 = integrate(PIEnsembleProblem(prob, pi2, (x, y) -> q₀ .+ ρ .* (2x - 1, 2y - 1)), DVRK(Gauss(2)))
nothing # hide
```

```@example lotka
grid = SecondPI{Float64, 2}(twoform, (15, 15), SecondFinDiffPlan)
solg = integrate(PIEnsembleProblem(prob, grid, (x, y) -> q₀ .+ ρ .* (2x - 1, 2y - 1)), DVRK(Gauss(2)))

plot_surface(grid, solg; xlabel = "q₁", ylabel = "q₂")
```

```@example lotka
plot_invariant(pi2, sol2; title = "Second Poincaré invariant")
```

The relative errors stay at the level of machine precision: with the correct form of the
symplectic potential, the `DVRK` variational integrator preserves the noncanonical Poincaré
invariants to machine accuracy, even though the Lotka-Volterra flow deforms the initial curve
and surface.
