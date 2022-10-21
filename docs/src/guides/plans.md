# Using Different Integral Implementations

There are a number of different implementations to choose from to calculate the invariants. For the first invariant, the choices are

```Julia
FirstPI{T, D}(θ, N, FirstFinDiffPlan)
FirstPI{T, D}(θ, N, FirstFourierPlan)
```

where `FirstFourierPlan` is the default if no plan type is given. `FirstFinDiffPlan` uses a finite difference approximation of the derivatives and calculates the final integral via the trapezoid rule. `FirstFourierPlan` transforms to frequency space, where the derivative is trivial, and computes the final integral in frequency space.

For the second invariant, there are two choices:

```Julia
SecondPI{T, D}(θ, N, SecondChebyshevPlan)
SecondPI{T, D}(θ, (Nx, Ny), SecondFinDiffPlan)
```

where `SecondChebyshevPlan` is the default. `SecondFinDiffPlan` uses finite differences to calculate the derivatives of the parameterisation and calculates the final integral via Simpson's rule. If the parameterisation is a second order polynomial, this method is exact. Note, that the finite difference implementation can use an arbitrary grid, which maybe specified with a tuple `(Nx, Ny)`. Alternatively, giving just a number `N` will result in using a square grid of points. `SecondChebyshevPlan` uses the Padua transform, as implemented in `ChebyshevTransforms.jl` to transform to a Chebyshev polynomial basis. In this basis, derivatives are calculated and then transformed back again. After evaluating the differential form at all points and contracting with the derivatives, another Padua transform is used to evaluate the final integral.
