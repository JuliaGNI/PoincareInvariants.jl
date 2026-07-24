# Harmonic Oscillator

The harmonic oscillator is the simplest canonical Hamiltonian system. With unit mass and
spring constant $k$ its Hamiltonian is

```math
H(q, p) = \frac{p^2}{2} + \frac{k \, q^2}{2} ,
```

so that the equations of motion in the canonical phase space $(q, p)$ read

```math
\dot{q} = p , \qquad \dot{p} = -k \, q .
```

We use the ready-made problem from
[GeometricProblems.jl](https://github.com/JuliaGNI/GeometricProblems.jl). For each invariant
we show the advected curve or surface, and two plots of the relative error over time: one for
the symplectic partitioned method `PartitionedGauss` (applied to the canonical `HODEProblem`),
and one comparing the symplectic `ImplicitMidpoint` method with the non-symplectic explicit
Runge-Kutta method `RK4` (both applied to the plain `ODEProblem` whose state is the pair
$[q, p]$).

```@example oscillator
using PoincareInvariants
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using CairoMakie

probh = hodeproblem([0.0], [0.0]; timespan = (0.0, 5.0), timestep = 0.2)  # canonical
probo = odeproblem([0.0, 0.0];   timespan = (0.0, 5.0), timestep = 0.2)   # plain ODE

saved_times(sol) = [sol[1].t[n] for n in 0:ntime(sol[1])]
relerr(Is)       = (Is .- Is[1]) ./ Is[1]
tsteps(sol)      = round.(Int, range(0, ntime(sol[1]), length = 10))  # 10 time slices

# (q, p) phase space points of every ensemble member at saved time index n
phasepoints(sol, n) = ([sol[j].q[n][1] for j in 1:nsamples(sol)],
                       [sol[j].p[n][1] for j in 1:nsamples(sol)])

# 3D view shared by all advected curve/surface plots
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

## First invariant

The first Poincaré invariant is the loop integral of the canonical one-form $\vartheta = p \, dq$,

```math
I_{1} = \oint_{\gamma} p \, dq ,
```

which equals the area enclosed by the closed curve $\gamma$. We take a circle of radius
$r = 1/2$ centred at the origin and integrate it with the three methods.

```@example oscillator
r = 0.5
init1 = ϕ -> (r * sinpi(2ϕ), r * cospi(2ϕ))

pi1 = CanonicalFirstPI{Float64, 2}(500)

sol1_pg = integrate(PIEnsembleProblem(probh, pi1, init1), PartitionedGauss(1))
sol1_im = integrate(PIEnsembleProblem(probo, pi1, init1), ImplicitMidpoint())
sol1_rk = integrate(PIEnsembleProblem(probo, pi1, init1), RK4())

Is1_pg = compute!(pi1, sol1_pg)
Is1_im = compute!(pi1, sol1_im)
Is1_rk = compute!(pi1, sol1_rk)

ts = saved_times(sol1_pg)
nothing # hide
```

As the loop is advected by the flow it rotates rigidly in phase space (the harmonic
oscillator is linear):

```@example oscillator
fig = Figure()
ax = Axis3(fig[1, 1]; xlabel = "q", ylabel = "p", zlabel = "t", view3..., title = "Advected Loop")
plot_loop!(ax, sol1_pg)
fig
```

The symplectic partitioned method conserves the invariant to machine precision:

```@example oscillator
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "Relative Error (I₁(t)-I₁(0))/I₁(0)",
    title = "PartitionedGauss")
hlines!(ax, [0.0]; color = :gray, linestyle = :dash)
scatter!(ax, ts, relerr(Is1_pg))
xlims!(ax, first(ts), last(ts))
fig
```

`ImplicitMidpoint` is symplectic too and likewise conserves the invariant to machine
precision, while the non-symplectic explicit Runge-Kutta method lets it drift:

```@example oscillator
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "Relative Error (I₁(t)-I₁(0))/I₁(0)",
    title = "ImplicitMidpoint vs Explicit Runge-Kutta-4")
hlines!(ax, [0.0]; color = :gray, linestyle = :dash)
scatter!(ax, ts, relerr(Is1_im); label = "ImplicitMidpoint")
scatter!(ax, ts, relerr(Is1_rk); label = "Explicit Runge-Kutta-4")
xlims!(ax, first(ts), last(ts))
axislegend(ax; position = :lb)
fig
```

## Second invariant

The second Poincaré invariant is the surface integral of the canonical two-form
$\omega = dq \wedge dp$,

```math
I_{2} = \int_{S} dq \wedge dp ,
```

i.e. the (signed) area of the surface $S$. We take a unit square centred at the origin.

```@example oscillator
init2 = (x, y) -> (x - 0.5, y - 0.5)

pi2 = CanonicalSecondPI{Float64, 2}(2_000)

sol2_pg = integrate(PIEnsembleProblem(probh, pi2, init2), PartitionedGauss(1))
sol2_im = integrate(PIEnsembleProblem(probo, pi2, init2), ImplicitMidpoint())
sol2_rk = integrate(PIEnsembleProblem(probo, pi2, init2), RK4())

Is2_pg = compute!(pi2, sol2_pg)
Is2_im = compute!(pi2, sol2_im)
Is2_rk = compute!(pi2, sol2_rk)
nothing # hide
```

For the surface we advect a coarse regular grid (with the same symplectic method) so it can
be drawn as a filled patch at each time:

```@example oscillator
grid = CanonicalSecondPI{Float64, 2}((15, 15), SecondFinDiffPlan)
nx, ny = getpointspec(grid)
solg = integrate(PIEnsembleProblem(probh, grid, init2), PartitionedGauss(1))

fig = Figure()
ax = Axis3(fig[1, 1]; xlabel = "q", ylabel = "p", zlabel = "t", view3..., title = "Advected Surface")
plot_surface!(ax, solg, nx, ny)
fig
```

```@example oscillator
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "Relative Error (I₂(t)-I₂(0))/I₂(0)",
    title = "PartitionedGauss")
hlines!(ax, [0.0]; color = :gray, linestyle = :dash)
scatter!(ax, ts, relerr(Is2_pg))
xlims!(ax, first(ts), last(ts))
fig
```

```@example oscillator
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "Relative Error (I₂(t)-I₂(0))/I₂(0)",
    title = "ImplicitMidpoint vs Explicit Runge-Kutta-4")
hlines!(ax, [0.0]; color = :gray, linestyle = :dash)
scatter!(ax, ts, relerr(Is2_im); label = "ImplicitMidpoint")
scatter!(ax, ts, relerr(Is2_rk); label = "Explicit Runge-Kutta-4")
xlims!(ax, first(ts), last(ts))
axislegend(ax; position = :lb)
fig
```

As for the first invariant, the two symplectic methods conserve the second invariant to
machine precision, whereas the explicit Runge-Kutta method does not.
