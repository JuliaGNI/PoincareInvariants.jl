# Poincaré Integral Invariant Example


## 2nd Poincaré Invariant

Load modules
```@example 2
using GeometricIntegrators
using PoincareInvariants
using PyPlot
```

Define functions for phasespace surface, one- and two-form and vector field
```@example 2
function surface(s, t, offset=1.0, factor=1.0)
    [(offset+s)/factor, (offset+t)/factor, s, t]
end

function theta(t, q)
    return [x[3] - x[2], x[4], 0., 0.]
end

function omega(t, q, B)
    B .=   [[ 0.  1. -1.  0.]
            [-1.  0.  0. -1.]
            [ 1.  0.  0.  0.]
            [ 0.  1.  0.  0.]]
end

function f(t, q, f)
    f[1] =  q[3]
    f[2] =  q[4]
    f[3] = +q[4]
    f[4] = -q[3]
end
nothing # hide
```

Set parameters
```@example 2
const nd = 4
const n₀ = 100
const nt = 200
const Δt = 0.1
nothing # hide
```


## 2nd Poincaré Invariant • Trapezoidal Quadrature

Create Poincaré Invariant
```@example 2
pinv = PoincareInvariant2ndTrapezoidal(
            x₀ -> ode = ODE(f, x₀), surface, omega,
            Δt, nd, n₀, n₀, nt)
nothing # hide
```

Compute solution
```@example 2
tab = TableauGLRK(1)
int = Integrator(pinv.equ, tab, pinv.Δt)
sol = integrate(pinv.equ, int, pinv.ntime)
nothing # hide
```

Plot solution
```@example 2
fig = figure(figsize=(4,4))
plot(sol.q[1,0,:], sol.q[2,0,:], ".", markersize=.5)
savefig("example_2nd_trapezoidal_initial_state.svg")

fig = figure(figsize=(4,4))
plot(sol.q[1,end,:], sol.q[2,end,:], ".", markersize=.5)
savefig("example_2nd_trapezoidal_final_state.svg")

fig = figure(figsize=(4,4))
for i in 1:n₀
    plot3D(sol.q[1,:,i], sol.q[2,:,i], collect(0:nt)*Δt)
end
savefig("example_2nd_trapezoidal_evolution1.png")

fig = figure(figsize=(4,4))
for i in 0:20:nt
    plot3D(sol.q[1,i,:], sol.q[2,i,:], i*Δt)
end
savefig("example_2nd_trapezoidal_evolution2.png")
```

![Initial State](example_2nd_trapezoidal_initial_state.svg)
![Final State](example_2nd_trapezoidal_final_state.svg)

![Evolution](example_2nd_trapezoidal_evolution1.png)
![Evolution](example_2nd_trapezoidal_evolution2.png)

Compute Poincarè invariant
```@example 2
I, J, ΔI, ΔJ = evaluate_poincare_invariant(pinv, sol)
nothing # hide
```

Plot invariant error
```@example 2
yf = matplotlib[:ticker][:ScalarFormatter]()
yf[:set_powerlimits]((-1,+1))
yf[:set_scientific](true)
yf[:set_useOffset](true)

fig = figure(figsize=(8,4))
plot((0:nt)*Δt, ΔI)
ax = gca()
ax[:yaxis][:set_major_formatter](yf)
savefig("example_2nd_trapezoidal.svg")
```

![](example_2nd_trapezoidal.svg)


## 2nd Poincaré Invariant • ApproxFun

Create Poincaré Invariant
```@example 2
pinv = PoincareInvariant2ndApproxFun(
            x₀ -> ode = ODE(f, x₀), surface, omega,
            Δt, nd, n₀, n₀, nt)
nothing # hide
```

Compute solution
```@example 2
tab = TableauGLRK(1)
int = Integrator(pinv.equ, tab, pinv.Δt)
sol = integrate(pinv.equ, int, pinv.ntime)
nothing # hide
```

Plot solution
```@example 2
fig = figure(figsize=(4,4))
plot(sol.q[1,0,:], sol.q[2,0,:], ".", markersize=.5)
savefig("example_2nd_approxfun_initial_state.svg")

fig = figure(figsize=(4,4))
plot(sol.q[1,end,:], sol.q[2,end,:], ".", markersize=.5)
savefig("example_2nd_approxfun_final_state.svg")

fig = figure(figsize=(4,4))
for i in 1:n₀
    plot3D(sol.q[1,:,i], sol.q[2,:,i], collect(0:nt)*Δt)
end
savefig("example_2nd_approxfun_evolution1.png")

fig = figure(figsize=(4,4))
for i in 0:20:nt
    plot3D(sol.q[1,i,:], sol.q[2,i,:], i*Δt)
end
savefig("example_2nd_approxfun_evolution2.png")
```

![Initial State](example_2nd_approxfun_initial_state.svg)
![Final State](example_2nd_approxfun_final_state.svg)

![Evolution](example_2nd_approxfun_evolution1.png)
![Evolution](example_2nd_approxfun_evolution2.png)


Compute Poincarè invariant
```@example 2
I, J, ΔI, ΔJ = evaluate_poincare_invariant(pinv, sol)
nothing # hide
```

Plot invariant error
```@example 2
yf = matplotlib[:ticker][:ScalarFormatter]()
yf[:set_powerlimits]((-1,+1))
yf[:set_scientific](true)
yf[:set_useOffset](true)

fig = figure(figsize=(8,4))
plot((0:nt)*Δt, ΔI)
ax = gca()
ax[:yaxis][:set_major_formatter](yf)
savefig("example_2nd_approxfun.svg")
```

![](example_2nd_approxfun.svg)


## 2nd Poincaré Invariant • OrthogonalPolynomialsQuasi


