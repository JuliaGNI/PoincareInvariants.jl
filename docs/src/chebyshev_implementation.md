# Chebyshev Implementation

Here, I describe how the second Poincaré invariant is calculated using the `ChebyshevPlan`.
The code for the `compute!` function in this implementation looks something like:

```Julia
function compute!(plan::ChebyshevPlan, Ω::Callable, phasepoints, t, p)
    paduatransform!(plan.phasecoeffs, plan.paduaplan, phasepoints)
    differentiate!(plan.∂x, plan.∂y, plan.diffplan, plan.phasecoeffs)
    getintegrand!(plan.intcoeffs, plan.intplan, Ω, phasepoints, t, p, plan.∂x, plan.∂y)
    integrate(plan.intcoeffs, plan.intweights)
end
```

First, we approximate the surface using Chebyshev polynomials. Coefficients are calculated use the Padua transform, which yields a matrix of coefficients. (See docs on PaduaTransforms)

Second, we differentiate the approximated function. Chebyshev polynomials can be differentiated via a linear transformation of the coefficients, yielding coefficients of a new polynomial in the Chebyshev basis. The linear transformation is applied to each column to differentiate with respect to one parameterisation variable and then to each row to differentiate with respect to the other. (See the functions `getdiffmat!` and `differentiate!` in the source code)

Third, the approximation of the derivatives is evaluated at the Padua points using the inverse Padua transform. The inital forward transform and the two inverse transforms together make up the vast majority of the computational cost.

Fourth, at each Padua point the derivatives in both directions (two vectors) and the invariant two-form are contracted via a vector matrix vector product. This is the integrand evaluated at each Padua point. (See `getintegrand!`)

Finally, we put the integrand values through one final Padua transform, obtaining a polynomial approximation of the integrand. This is then integrated in 2D by another vector matrix vector product, where the matrix is the matrix of Chebyshev coefficients and the vector is a vector of definite integrals of each of the Chebyshev polynomials. (See `getintweights` for the calculation of the definite integrals and `integrate` to perform the integration.)
