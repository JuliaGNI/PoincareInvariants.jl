# Implementation

Evaluating the invariants numerically (see [Theory](theory.md)) requires sampling the advected loop or surface, differentiating the parametrisation, evaluating the differential form, and integrating the result. To obtain high accuracy — which is essential for chaotic systems, where the sampling of the submanifold degrades quickly — the default plans use *spectral* methods, whose derivatives and quadratures converge exponentially for smooth data. The geometry of the two submanifolds dictates the choice of basis: the first invariant lives on a periodic loop and is handled with the fast Fourier transform, while the second invariant lives on a square and is handled with the Padua transform to a Chebyshev basis, as implemented in [`ChebyshevTransforms.jl`](https://juliagni.github.io/ChebyshevTransforms.jl/).

## First invariant via the Fourier transform

The first invariant is computed by `FirstFourierPlan` in `src/FirstFourierPlans.jl`. Recall
```math
I_1 (t) = \int_0^1 \vartheta_i (q_{(\tau)}(t)) \, \frac{d q_{(\tau)}^i}{d\tau} \, d\tau ,
```
where the loop is closed and hence periodic in the parameter ``\tau``. The loop is sampled at ``N`` equidistant points ``\tau_j = j/N``, ``j = 0, \dots, N-1`` (`getpoints`, `src/FirstFourierPlans.jl:70`), which is the natural grid for a periodic domain: on such a domain the trapezoidal rule is spectrally accurate ([Trefethen & Weideman, 2014](https://doi.org/10.1137/130932132)).

Because the data are periodic, both the coordinate samples ``q^i(\tau_j)`` and the one-form component samples ``\vartheta_i(q(\tau_j))`` are expanded in a discrete Fourier series. In frequency space differentiation is a diagonal operation — the ``k``-th Fourier coefficient of ``dq^i/d\tau`` is simply ``2\pi i k`` times the coefficient of ``q^i`` — so no finite differencing is needed. The integral itself is then evaluated in frequency space via Parseval's relation, avoiding any explicit quadrature over the physical grid.

Concretely, `compute!` (`src/FirstFourierPlans.jl:33`) evaluates the one-form at every sample, and for each phase-space dimension ``d`` applies a real FFT (an `FFTW.rFFTWPlan`) to both the one-form samples and the coordinate samples, accumulating
```math
I_1 \approx \sum_{d} \sum_{k} \operatorname{Re} \left[ \hat{z}_k \, \overline{\hat{\vartheta}_k} \; \frac{4 \pi i \, (k-1)}{N^2} \right] ,
```
which corresponds to `src/FirstFourierPlans.jl:55-57`. Here ``\hat{z}_k`` and ``\hat{\vartheta}_k`` are the transformed coordinate and one-form coefficients, the factor ``i\,(k-1)`` implements the derivative (with ``k-1`` the integer wavenumber under Julia's 1-based indexing), and the factor ``4\pi/N^2`` combines the normalisation of the unnormalised FFT with the factor of two arising from the Hermitian symmetry of a real signal, whose real FFT stores only the non-negative frequencies ``k = 1, \dots, N/2 + 1``.

## Second invariant via the Padua transform

The second invariant is computed by `SecondChebyshevPlan` in `src/SecondChebyshevPlans.jl`, which represents the surface by Chebyshev polynomials. Recall
```math
I_2 (t) = \int_0^1 \int_0^1 \omega_{ij} (q_{(\sigma, \tau)}(t)) \, \frac{d q_{(\sigma, \tau)}^i}{d\sigma} \, \frac{d q_{(\sigma, \tau)}^j}{d\tau} \, d\sigma \, d\tau .
```
A Chebyshev representation allows an extremely accurate approximation of the surface even when it becomes severely deformed. The computation proceeds in four stages, orchestrated by `compute!` (`src/SecondChebyshevPlans.jl:155`).

### Sampling at the Padua points

The surface is sampled at the *Padua points*, the (near-)optimal nodal set for bivariate polynomial interpolation on the square. `getpoints` (`src/SecondChebyshevPlans.jl:177`) calls `getpaduapoints` from `ChebyshevTransforms.jl`, mapping the parameter domain ``[0,1]^2`` onto the Chebyshev domain ``[-1,1]^2``. The requested point count ``N`` is rounded up to a valid Padua number via `getpointspec` / `nextpaduanum` (`src/SecondChebyshevPlans.jl:175`), which fixes the total polynomial degree.

### Transform to the Chebyshev basis

`paduatransform!` transforms the sampled coordinates into Chebyshev coefficients ``a_{ij}`` of the expansion
```math
q(\sigma, \tau) \approx \sum_{i, j} a_{ij} \, T_i(\tau) \, T_j(\sigma) ,
```
where ``T_n(x) = \cos(n \arccos x)`` is the ``n``-th Chebyshev polynomial. As detailed in the [Padua Transforms documentation](https://juliagni.github.io/ChebyshevTransforms.jl/stable/padua_transforms/), the transform is implemented as a discrete cosine transform (itself an FFT), which is what makes it fast.

### Differentiation in the Chebyshev basis

Derivatives are cheap in the Chebyshev basis: they amount to multiplication by a fixed triangular differentiation matrix (`getdiffmat`, `src/SecondChebyshevPlans.jl:19`). `differentiate!` (`src/SecondChebyshevPlans.jl:42`) obtains ``\partial_\sigma`` by right-multiplying the coefficient matrix and ``\partial_\tau`` by left-multiplying it with the (transposed) differentiation matrix.

### Contraction and Clenshaw–Curtis integration

The derivative coefficients are transformed back to values at the Padua points with `invpaduatransform!`, the two-form ``\omega`` is evaluated at each point, and the integrand ``\omega_{ij}\,(\partial_\sigma q^i)(\partial_\tau q^j)`` is assembled pointwise as `dot(∂y, ω, ∂x)` (`getintegrand!`, `src/SecondChebyshevPlans.jl:95-119`). A further `paduatransform!` expresses this integrand in the Chebyshev basis, whose coefficients are then integrated exactly using the Chebyshev integration weights
```math
w_i = \int_{-1}^{1} T_i(x) \, dx = \begin{cases} \dfrac{2}{1 - i^2} & i \ \text{even} \\[1ex] 0 & i \ \text{odd} \end{cases}
```
(`getintweights`, `src/SecondChebyshevPlans.jl:68`). The final value is the double contraction ``\sum_{ij} w_i \, a_{ij} \, w_j``, evaluated as `integrate(coeffs, intweights) = dot(intweights, coeffs, intweights)` (`src/SecondChebyshevPlans.jl:71`) — a two-dimensional Clenshaw–Curtis cubature.

Both plans have finite-difference counterparts (`FirstFinDiffPlan`, `SecondFinDiffPlan`) for cases where a spectral basis is unsuitable; see [Using Different Integral Implementations](guides/plans.md).
