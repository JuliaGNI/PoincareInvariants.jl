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
Poincaré invariants — to machine accuracy.

```@example lotka
using PoincareInvariants
using GeometricIntegrators
using GeometricProblems.LotkaVolterra2dSingular
using StaticArrays
using CairoMakie

prob = lodeproblem([2.0, 1.0]; timespan = (0.0, 1.0), timestep = 0.02)

saved_times(sol) = [sol[1].t[n] for n in 0:ntime(sol[1])]
relerr(Is)       = (Is .- Is[1]) ./ Is[1]
tsteps(sol)      = round.(Int, range(0, ntime(sol[1]), length = 10))

phasepoints(sol, n) = ([sol[j].q[n][1] for j in 1:nsamples(sol)],
                       [sol[j].q[n][2] for j in 1:nsamples(sol)])

view3   = (; azimuth = 1.775π, elevation = π / 8, aspect = (1, 1, 1.6))
palette = Makie.wong_colors()

function plot_loop!(ax, sol)
    ts = saved_times(sol)
    for (k, n) in enumerate(tsteps(sol))
        xs, ys = phasepoints(sol, n)
        lines!(ax, xs, ys, fill(ts[n + 1], length(xs)); color = palette[mod1(k, 7)])
    end
    zlims!(ax, first(ts), last(ts))
end

function plot_surface!(ax, sol, nx, ny)
    ts = saved_times(sol)
    for (k, n) in enumerate(tsteps(sol))
        xs, ys = phasepoints(sol, n)
        surface!(ax, reshape(xs, ny, nx), reshape(ys, ny, nx), fill(ts[n + 1], ny, nx);
            color = fill(palette[mod1(k, 7)], ny, nx), shading = NoShading)
    end
    zlims!(ax, first(ts), last(ts))
end
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
Is1  = compute!(pi1, sol1)

ts = saved_times(sol1)
nothing # hide
```

The loop is transported and deformed by the Lotka-Volterra flow:

```@example lotka
fig = Figure()
ax = Axis3(fig[1, 1]; xlabel = "q₁", ylabel = "q₂", zlabel = "t", view3..., title = "Advected Loop")
plot_loop!(ax, sol1)
fig
```

```@example lotka
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "Relative Error (I₁(t)-I₁(0))/I₁(0)",
    title = "First Poincaré invariant")
hlines!(ax, [0.0]; color = :gray, linestyle = :dash)
scatter!(ax, ts, relerr(Is1))
xlims!(ax, first(ts), last(ts))
fig
```

## Second invariant

```math
I_{2} = \int_{S} \omega_{ij}(q) \, dq^i \, dq^j
```

We advect a square around $(2, 1)$, staying within the positive quadrant.

```@example lotka
pi2  = SecondPI{Float64, 2}(twoform, 2_000)
sol2 = integrate(PIEnsembleProblem(prob, pi2, (x, y) -> q₀ .+ ρ .* (2x - 1, 2y - 1)), DVRK(Gauss(2)))
Is2  = compute!(pi2, sol2)
nothing # hide
```

```@example lotka
grid = SecondPI{Float64, 2}(twoform, (15, 15), SecondFinDiffPlan)
nx, ny = getpointspec(grid)
solg = integrate(PIEnsembleProblem(prob, grid, (x, y) -> q₀ .+ ρ .* (2x - 1, 2y - 1)), DVRK(Gauss(2)))

fig = Figure()
ax = Axis3(fig[1, 1]; xlabel = "q₁", ylabel = "q₂", zlabel = "t", view3..., title = "Advected Surface")
plot_surface!(ax, solg, nx, ny)
fig
```

```@example lotka
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "Relative Error (I₂(t)-I₂(0))/I₂(0)",
    title = "Second Poincaré invariant")
hlines!(ax, [0.0]; color = :gray, linestyle = :dash)
scatter!(ax, ts, relerr(Is2))
xlims!(ax, first(ts), last(ts))
fig
```

The relative errors stay at the level of machine precision: with the correct form of the
symplectic potential, the `DVRK` variational integrator preserves the noncanonical Poincaré
invariants to machine accuracy, even though the Lotka-Volterra flow deforms the initial curve
and surface.
