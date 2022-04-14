# Padua Transforms

```@meta
CurrentModule = PoincareInvariants.SecondChebyshevPlans.PaduaTransforms
DocTestSetup = quote
    using PoincareInvariants.SecondChebyshevPlans.PaduaTransforms
end
```

## Introduction

The Padua transform yields coefficients ``a_{ij}`` used to approximate a bivariate function ``f(x, y)`` with

```math
f(x, y) ≈ \sum_{i, j} a_{ij} \; T_i(y) \; T_j(x)
```

where ``T_n(x) := \cos(n \arccos(x))`` is the nth Chebyshev polynomial. For reference, the first few Chebyshev polynomials are

```math
\begin{aligned}
T_0(x) &= 1 \\
T_1(x) &= x \\
T_2(x) &= 2x^2 - 1 \\
T_3(x) &= 4x^3 - 3x \\
T_4(x) &= 8x^4 - 8x^2 + 1
\end{aligned}
```

The inverse transform takes the coefficients ``a_{ij}`` and returns the function ``f(x, y)`` evaluated at the so called Padua points.

## Basic Usage

Start by evaluating a function on the Padua points. Here we approximate the function with a polynomial of total degree 3.
The total degree is the degree of the largest polynomial in x plus the degree of the largest polynomial in y.

```jldoctest basic
julia> vals = getpaduapoints(3) do x, y
           y * (2x^2 - 1) + 5 * x * y + 2.5
       end
10-element Vector{Float64}:
  8.5
  2.5
 -3.5
  3.914213562373095
  1.0857864376269049
 -0.5
  2.5
  5.5
 -0.3284271247461903
  5.32842712474619
```

Then create a [`PaduaTransformPlan`](@ref) and apply the [`paduatransform!`](@ref) to get the Chebyshev coefficeints ``a_{ij}``.

```jldoctest basic
julia> plan = PaduaTransformPlan{Float64}(3);

julia> coeffs = paduatransform!(zeros(4, 4), plan, vals)
4×4 Matrix{Float64}:
 2.5  0.0  0.0  0.0
 0.0  5.0  1.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
```

We can go back to values using an [`InvPaduaTransformPlan`](@ref) and applying [`invpaduatransform!`](@ref)

```jldoctest basic
julia> invplan = InvPaduaTransformPlan{Float64}(3);

julia> out = invpaduatransform!(Vector{Float64}(undef, getpaduanum(3)), invplan, coeffs)
10-element Vector{Float64}:
  8.5
  2.5
 -3.5
  3.914213562373095
  1.0857864376269049
 -0.5
  2.5
  5.5
 -0.3284271247461903
  5.32842712474619

julia> out ≈ vals
true
```

where we used [`getpaduanum`](@ref) to get the number of Padua points and coefficients corresponding to total degree 3.

## The Algorithm -- Step by Step

To obtain the coefficients ``a_{ij}`` we need to evaluate the function ``f`` at some points ``(x_l, y_k)`` and evaluate

```math
a_{ij} = \sum_{k, l} \; f(x_l, y_k) \; T_i(y_k) \; T_j(x_l)
```

If we let ``x_l = \cos{\frac{lπ}{n}}`` and ``y_k = \cos{\frac{kπ}{m}}``, we have

```math
a_{ij} = \sum_{k, l} \; f(\cos{\frac{lπ}{n}}, \cos{\frac{kπ}{m}}) \; T_i(\cos{\frac{kπ}{m}}) \; T_j(\cos{\frac{lπ}{n}})
```

which simplifies to

```math
a_{ij} = \sum_{k, l} \; f(\cos{\frac{lπ}{n}}, \cos{\frac{kπ}{m}}) \; \cos{\frac{ikπ}{m}} \; \cos{\frac{jlπ}{n}}
```

because of the definition of the Chebyshev polynomials as ``T_n(x) := \cos(n \arccos(x))``.
The expression above looks like a discrete cosine transform, which is what we will use to implement the Padua transform.
Note, that the formula above is not quite exact since we need to apply a weighting factor to the coefficients.
For further reading on the Padua Transform, please see the following:

For details on the Padua points as good nodes for polynomial interpolation:
[Marco Caliari, Stefano De Marchi, Marco Vianello. Bivariate polynomial interpolation on the square at new nodal sets](https://doi.org/10.1016/j.amc.2004.07.001)

and for details on the implementation of the transform via a discrete cosine transform:
[Marco Caliari, Stefano De Marchi, Alvise Sommariva, Marco Vianello. Padua2DM: fast interpolation and cubature at the Padua points in Matlab/Octave](https://link.springer.com/content/pdf/10.1007/s11075-010-9373-1.pdf)

### The Padua Points

There are multiple ways to define the set of Padua points. For our purposes we will use

```math
\textrm{Pad}_n = \{(\cos{\frac{jπ}{n}}, \cos{\frac{iπ}{n + 1}}) \; | \; 0 ≤ j ≤ n, \; 0 ≤ i ≤ n + 1, \; i - j \; \textrm{even} \}
```

where ``n`` is the total degree of the polynomial we can approximate using the Padua points.
To generate the Padua points we can use the function [`getpaduapoints`](@ref) as follows

```jldoctest
julia> getpaduapoints(3)
10×2 Matrix{Float64}:
  1.0   1.0
  1.0   0.0
  1.0  -1.0
  0.5   0.707107
  0.5  -0.707107
 -0.5   1.0
 -0.5   0.0
 -0.5  -1.0
 -1.0   0.707107
 -1.0  -0.707107
```

Using do block syntax we can evaluate a function on the Padua points.

```jldoctest step_by_step
julia> vals = getpaduapoints(3) do x, y
                  y * (2x^2 - 1) + 5 * x * y + 2.5
              end
10-element Vector{Float64}:
  8.5
  2.5
 -3.5
  3.914213562373095
  1.0857864376269049
 -0.5
  2.5
  5.5
 -0.3284271247461903
  5.32842712474619
```

### The Padua Transform

Having evaluated our function on the Padua points, we can perform the transform. First, we initialise a transform plan.

```jldoctest step_by_step
julia> plan = PaduaTransformPlan{Float64}(3);
```

The transform consists of writing `vals` into `plan.vals`, applying a fast fourier transform and weighting the resulting coefficients to obtain the Chebyshev coefficients. We start by writting values into the values matrix using [`tovalsmat!`](@ref).

```jldoctest step_by_step
julia> PaduaTransforms.tovalsmat!(plan.vals, vals, 3)
5×4 Matrix{Float64}:
  8.5  0.0      -0.5   0.0
  0.0  3.91421   0.0  -0.328427
  2.5  0.0       2.5   0.0
  0.0  1.08579   0.0   5.32843
 -3.5  0.0       5.5   0.0
```

`vals` are written into the matrix `plan.vals` such that the entry `plan.vals[i+1, j+1]` corresponds
to the Padua point ``(\cos{\frac{jπ}{n}}, \cos{\frac{iπ}{n + 1}})``. However, since the Padua points are only those points with
``i-j`` even, all entries corresponding to ``i-j`` odd are left out. These off grid entries must be filled with 0.
Next, we can apply the discrete cosine transform and [`weight!`](@ref) the coefficients to obtain the Chebyshev coefficients.

```jldoctest step_by_step
julia> plan.dctplan * plan.vals
5×4 Matrix{Float64}:
 60.0   0.0   0.0   0.0
  0.0  30.0   6.0   0.0
  0.0   0.0   0.0   0.0
  0.0   6.0  30.0   0.0
  0.0   0.0   0.0  60.0

julia> PaduaTransforms.weight!(plan.vals, 3)
5×4 Matrix{Float64}:
 2.5  0.0  0.0  0.0
 0.0  5.0  1.0  0.0
 0.0  0.0  0.0  0.0
 0.0  1.0  5.0  0.0
 0.0  0.0  0.0  2.5
```

The weighting factor we apply to the coefficients is

```math
w = \frac{1}{n(n+1)} ⋅ \begin{cases}
    \frac{1}{2} & \textrm{if on vertex}   \\
    1           & \textrm{if on edge}     \\
    2           & \textrm{if in interior} \\
\end{cases}
```

Finally, we write those coefficients corresponding to a total degree of 3 or lower (the upper left triangular) into the output.

```jldoctest step_by_step
julia> coeffs = zeros(4, 4)
4×4 Matrix{Float64}:
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0

julia> PaduaTransforms.fromcoeffsmat!(coeffs, plan.vals, 3)
4×4 Matrix{Float64}:
 2.5  0.0  0.0  0.0
 0.0  5.0  1.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
```

### The Inverse Padua Transform

If we have a set of Chebyshev coefficients and want to obtain the values of the corresponding Chebyshev polynomial on the Padua points,
we must use the inverse Padua transform. Again, we start with a plan and write our coefficients into `plan.coeffs`.

```jldoctest step_by_step
julia> invplan = InvPaduaTransformPlan{Float64}(3);

julia> PaduaTransforms.tocoeffsmat!(invplan.coeffs, coeffs)
5×4 Matrix{Float64}:
 2.5  0.0  0.0  0.0
 0.0  5.0  1.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
```

Next, we weight the coefficients using [`invweight!`](@ref) and apply a discrete cosine transform.

```jldoctest step_by_step
julia> PaduaTransforms.invweight!(invplan.coeffs)
5×4 Matrix{Float64}:
 2.5  0.0   0.0   0.0
 0.0  1.25  0.25  0.0
 0.0  0.0   0.0   0.0
 0.0  0.0   0.0   0.0
 0.0  0.0   0.0   0.0

julia> invplan.dctplan * invplan.coeffs
5×4 Matrix{Float64}:
  8.5      4.5      -0.5      -1.5
  6.74264  3.91421   0.37868  -0.328427
  2.5      2.5       2.5       2.5
 -1.74264  1.08579   4.62132   5.32843
 -3.5      0.5       5.5       6.5
```

[`invweight!`](@ref) applies the weighting

```math
w = \begin{cases}
    1           & \textrm{if on vertex}   \\
    \frac{1}{2} & \textrm{if on edge}     \\
    \frac{1}{4} & \textrm{if in interior} \\
\end{cases}
```

Finally, we copy over those values corresponding to ``i-j`` even and we have our values back.

```jldoctest step_by_step
julia> out = PaduaTransforms.fromvalsmat!(Vector{Float64}(undef, getpaduanum(3)), invplan.coeffs, 3)
10-element Vector{Float64}:
  8.5
  2.5
 -3.5
  3.914213562373095
  1.0857864376269049
 -0.5
  2.5
  5.5
 -0.3284271247461903
  5.32842712474619

julia> vals ≈ out
true
```

## Reference

```@docs
PaduaTransforms
```

### Padua Points

```@docs
getpaduapoints
```

### Numbers of Coefficients and Points

```@docs
getpaduanum
getdegree
nextpaduanum
nextdegree
```

### Padua Transform

```@docs
PaduaTransformPlan
paduatransform!
```

### Inverse Padua Transform

```@docs
InvPaduaTransformPlan
invpaduatransform!
```

### Internals

```@docs
paduapoint
ispadua
tovalsmat!
weight!
fromcoeffsmat!
tocoeffsmat!
invweight!
fromvalsmat!
```
