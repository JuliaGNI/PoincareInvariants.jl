# Massless Charged Particle

The massless charged particle in the plane is a noncanonical symplectic system. Its phase
space is just the position $x = (x_1, x_2)$, and its dynamics derive from the degenerate
Lagrangian

```math
L(x, \dot{x}) = A(x) \cdot \dot{x} - \phi(x) ,
```

with magnetic vector potential $A$, electrostatic potential $\phi$ and magnetic field
$B(x) = \nabla \times A(x)$. The Hamiltonian form of the equations of motion,

```math
\dot{x} = \frac{1}{B(x)} \begin{pmatrix} 0 & -1 \\ +1 & 0 \end{pmatrix} \nabla \phi(x) ,
```

carries a state-dependent (noncanonical) symplectic structure. We use the problem and its
exported forms from
[GeometricProblems.jl](https://github.com/JuliaGNI/GeometricProblems.jl), integrating the
degenerate `IODEProblem` with a **Degenerate Variational Runge-Kutta** method, `DVRK`, which
preserves the noncanonical structure.

```@example particle
using PoincareInvariants
using GeometricIntegrators
import GeometricProblems.MasslessChargedParticle as MCP
using StaticArrays
using CairoMakie

par  = MCP.default_parameters()
prob = MCP.iodeproblem([1.0, 1.0]; timespan = (0.0, 5.0), timestep = 0.1)

saved_times(sol) = [sol[1].t[n] for n in 0:ntime(sol[1])]
relerr(Is)       = (Is .- Is[1]) ./ Is[1]
tsteps(sol)      = round.(Int, range(0, ntime(sol[1]), length = 10))

# (x₁, x₂) positions of every ensemble member at saved time index n
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

The first invariant uses the one-form $\vartheta = A(x)$ (the magnetic vector potential), and
the second uses the two-form $\omega = d\vartheta$, which for this system is

```math
\omega = \begin{pmatrix} 0 & B(x) \\ -B(x) & 0 \end{pmatrix} .
```

Both are provided directly by the problem module (`ϑ` and `B`), so we only reshape them into
the `form(z, t, p)` interface expected by `PoincareInvariants`. The parameter slot `p`
carries the physical parameters `par`.

```@example particle
oneform(z, t, p) = SVector{2}(MCP.ϑ(t, z, p))

function twoform(z, t, p)
    b = MCP.B(z, p)
    @SMatrix [0.0  b;
              -b   0.0]
end
nothing # hide
```

## First invariant

```math
I_{1} = \oint_{\gamma} A(x) \cdot dx
```

We advect a small circle of radius $\rho = 0.2$ around the initial position $(1, 1)$ and pass
the parameters `par` to [`compute!`](@ref) so the form can evaluate $A$.

```@example particle
q₀ = SVector(1.0, 1.0)
ρ  = 0.2

pi1  = FirstPI{Float64, 2}(oneform, 500)
sol1 = integrate(PIEnsembleProblem(prob, pi1, ϕ -> q₀ .+ ρ .* (cospi(2ϕ), sinpi(2ϕ))), DVRK(Gauss(2)))
Is1  = compute!(pi1, sol1, par)

ts = saved_times(sol1)
nothing # hide
```

The loop is transported and deformed by the flow:

```@example particle
fig = Figure()
ax = Axis3(fig[1, 1]; xlabel = "x₁", ylabel = "x₂", zlabel = "t", view3..., title = "Advected Loop")
plot_loop!(ax, sol1)
fig
```

```@example particle
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
I_{2} = \int_{S} \omega = \int_{S} B(x) \, dx_1 \, dx_2
```

is the magnetic flux through the surface $S$. We advect a small square around $(1, 1)$.

```@example particle
pi2  = SecondPI{Float64, 2}(twoform, 2_000)
sol2 = integrate(PIEnsembleProblem(prob, pi2, (x, y) -> q₀ .+ ρ .* (2x - 1, 2y - 1)), DVRK(Gauss(2)))
Is2  = compute!(pi2, sol2, par)
nothing # hide
```

```@example particle
grid = SecondPI{Float64, 2}(twoform, (15, 15), SecondFinDiffPlan)
nx, ny = getpointspec(grid)
solg = integrate(PIEnsembleProblem(prob, grid, (x, y) -> q₀ .+ ρ .* (2x - 1, 2y - 1)), DVRK(Gauss(2)))

fig = Figure()
ax = Axis3(fig[1, 1]; xlabel = "x₁", ylabel = "x₂", zlabel = "t", view3..., title = "Advected Surface")
plot_surface!(ax, solg, nx, ny)
fig
```

```@example particle
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "Relative Error (I₂(t)-I₂(0))/I₂(0)",
    title = "Second Poincaré invariant")
hlines!(ax, [0.0]; color = :gray, linestyle = :dash)
scatter!(ax, ts, relerr(Is2))
xlims!(ax, first(ts), last(ts))
fig
```

The `DVRK(Gauss(2))` method is a fourth-order variational integrator for noncanonical
symplectic systems, so both invariants are conserved with high accuracy as the curve and
surface are advected by the flow. A standard, structure-agnostic integrator would let them
drift.
