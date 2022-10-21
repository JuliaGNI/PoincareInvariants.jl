# Compute Invariants With DifferentialEquations.jl

Before reading this, it's worth familiarising yourself with the [Ensemble Simulations](https://diffeq.sciml.ai/stable/features/ensemble/#ensemble) section of the DifferentialEquations.jl docs.

To compute invariants with DifferentialEquations.jl, we'll need to create a `PIEnsembleProblem`. This works as follows

```Julia
PIEnsembleProblem(init, prob, pinv;
    prob_func = (prob,i,repeat)->prob,
    output_func = (sol,i) -> (sol,false),
    reduction = (u,data,I)->(append!(u,data),false),
    u_init = [], safetycopy = false
)
```

`init` is the phasespace curve or surface parameterisation, which determines the initial conditions;
`prob` is a problem from DifferentialEquations.jl, which specifies the time span and differential equation;
`pinv` is an AbstractPoincareInvariant;
`prob_func` allows the user to `remake` the problem for each trajectory in the ensemble;
all further arguments are exactly like the `EnsembleProblem` type from DifferentialEquations.jl. You need not add a `prob_func` of your own to have the correct initial conditions. These will be entered automaticaly according ot the given setup objet and parameterisation. However, by default, the type of the initial condition is `Vector{T}`. If you want to use an `ArrayPartion` or a `StaticArray` or other type to represent a point in phase space, you should specify a `prob_func` to do the conversion.

All in all, this can look like so:

```Julia
using OrdinaryDiffEq
using RecursiveArrayTools: ArrayPartition
using PoincareInvariants

dt = 0.1
prob = SecondOrderODEProblem((p, θ, params, t) -> [-sin(θ[1])], 0.0, 0.0, (0.0, 2.0))
pf(prob, i, repeat) = remake(prob; u0 = ArrayPartition((prob.u0[1:1], prob.u0[2:2])))

pi1 = CanonicalFirstPI{Float64, 2}(1_000)
ens_prob = PIEnsembleProblem(ϕ -> (sinpi(2ϕ), 3 * cospi(2ϕ)), prob, pi1; prob_func=pf)
```

Once, we have our `PIEnsembleProblem`, we can call solve on it, just like the `EnsembleProblem` it wraps. The general pattern is

```Julia
solve(prob::PIEnsembleProblem, alg, ensemblealg; kwargs...)
```

Most arguments are simply passed to the solve function called on the wrapped `EnsembleProblem`.
The `trajectories` keyword argument is already set as `getpointnum(pinv)` and the `adaptive` keyword is set to false by default. Otherwise the behaviour is identical. For our example, we would have

```Julia
sol = solve(ens_prob, SymplecticEuler(), EnsembleSerial(); dt=dt)
```

Finally, we can call `compute!` on the result, which returns a `Vector{T}` of invariant values for every saved timestep. You must make sure all trajectories save at the same time steps.

```Julia
compute!(pi1, sol[, p])
```

The times and points are taken from the solution. An optional parameter `p` may be given, which is passed to the differential form.
