# Second Poincaré Invariants

```@meta
CurrentModule = PoincareInvariants.SecondPoincareInvariants
DocTestSetup = quote
    using PoincareInvariants
end
```

For this example, we will calculate the Second Poincaré Invariant for a simple example. We will be working in 2D phase space
for simplicity. To approximate our surface, we will use approximately 500 points.

```jldoctest usage
julia> D = 2; N = 500;
```

To calculate the second poincare invariant, we will need to define the two-form ``\Omega``, the integral invariant.

```jldoctest usage
julia> Ω(z, t, p) = CanonicalSymplecticMatrix(D);
```

The two-form is just a `D x D` matrix. In this case it is the canonical symplectic matrix.
Here, we have defined it as a function of the phase space position `z`, time `t` and arbitrary parameter `p`.
In future it will also be possible to use a constant matrix or an in-place function.

Now, we can initialise the setup object used to calculate the invariant.

```jldoctest usage
julia> pinv = SecondPoincareInvariant{Float64}(Ω, D, N);
```

The type `Float64` specifies that the final result as well as all intermediate calculations will use the type `Float64`.
The number `N` specifies the approximate number of points to be used to approximate the surface.
The exact number depends on implementation details. The type of `Ω` signals it is not a constant matrix.

The setup object contains the two-form, the dimension of the phase space and the number of points used to approximate the surface.
These properties can be probed using the functions [`getform`](@ref), [`getdim`](@ref) and [`getpointnum`](@ref).
The setup object also contains a `plan` object which is used to preallocate memory and setup the computation.

```jldoctest usage
julia> getform(pinv)
Ω (generic function with 1 method)

julia> getdim(pinv)
2

julia> getpointnum(pinv)
528
```

Let us use a half circle in phase space as our surface.
To get the points in phase space used in approximating our surface, we use the [`getpoints`](@ref) function.

```jldoctest usage
julia> phasepoints = getpoints(pinv) do r, θ
           s, c = sincospi(θ)
           return r .* (c, s)
       end;
```

Finally, we may compute the integral invariant using the function [`compute!`](@ref).
The last two arguments represent the time and any additional parameters. Both are passed to the two-form function `Ω(z, t, p)`.

```jldoctest usage
julia> p = compute!(pinv, phasepoints, 0, nothing)
1.570796326794897

julia> (p - π/2) < 5eps()
true
```

We can see that the answer is accurate within five times machine epsilon.
If we now evolve each point forward in time, we should see that the invariant is conserved.

```jldoctest usage
julia> function free_particle!(points, t)
           mid = length(points) ÷ 2
           for i in 1:mid
               points[i] .+= points[mid+i] .* t
           end
       end
free_particle! (generic function with 1 method)

julia> free_particle!(phasepoints, 100);

julia> abs(compute!(pinv, phasepoints, 0, nothing) - π/2) < 50eps()
true
```

We find that up to some numerical error the invariant is indeed conserved.
