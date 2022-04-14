"""
    PaduaTransforms

an implementation of the Padua transform and its inverse via the fast Fourier transform.
"""
module PaduaTransforms

using FFTW
using LinearAlgebra: rmul!

export getpaduanum, getdegree, nextdegree, nextpaduanum
export getpaduapoints
export PaduaTransformPlan, paduatransform!
export InvPaduaTransformPlan, invpaduatransform!

## Number of Padua Points and Degree ##

"""
    getpaduanum(n)

calculates number of Padua points needed to approximate a function using Chebyshev
polynomials up to total degree `n`. This number is equal to the number of coefficients.
The formula is

```math
N = (n + 1) ⋅ (n + 2) ÷ 2
```

# Examples
```jldoctest
julia> getpaduanum(13)
105
```
"""
getpaduanum(degree) = (degree + 1) * (degree + 2) ÷ 2

"""
    getdegree(N)

calculates total degree, given the number of coefficients or Padua points `N`.
Throws an error if `N` is not a possible number of Padua points. The formula is

```math
n = \\frac{\\sqrt{1 + 8N} - 3}{2}
```

# Examples
```jldoctest
julia> getdegree(105)
13

julia> getdegree(104)
ERROR: ArgumentError: number of Padua points or coeffs must be (n + 1) * (n + 2) ÷ 2
[...]
```
"""
function getdegree(N)
    d = (sqrt(1 + 8N) - 3) / 2
    isinteger(d) ? Int(d) : throw(ArgumentError(
        "number of Padua points or coeffs must be (n + 1) * (n + 2) ÷ 2"))
end

"""
    nextdegree(N)

get degree of `nextpaduanum(N)`.

# Examples

```jldoctest
julia> nextdegree(104)
13
```
"""
nextdegree(N) = Int(cld(sqrt(1 + 8N) - 3, 2))

"""
    nextpaduanum(N)

get next valid number of Padua points ≥ `N`.

# Examples

```jldoctest
julia> nextpaduanum(104)
105
```
"""
nextpaduanum(N) = getpaduanum(nextdegree(N))

## Padua Points ##

"""
    paduapoint(T::Type, j::Integer, i::Integer, n::Integer)

returns the Padua point ``z_{ij}``, where

```math
z_{ij} = (\\cos{\\frac{jπ}{n}}, \\cos{\\frac{iπ}{n+1}})
```

Note, that only points with ``i-j`` even are actually Padua points.
Check with [`ispadua`](@ref).

# Examples
```jldoctest
julia> [PaduaTransforms.paduapoint(Float32, x, y, 1) for y in 0:1+1, x in 0:1]
3×2 Matrix{Tuple{Float32, Float32}}:
 (1.0, 1.0)   (-1.0, 1.0)
 (1.0, 0.0)   (-1.0, 0.0)
 (1.0, -1.0)  (-1.0, -1.0)
```

"""
function paduapoint(::Type{T}, j, i, n) where T
    x = cospi(T(j) / T(n))
    y = cospi(T(i) / T(n + 1))
    return x, y
end

"""
    ispadua(i, j)

returns if [`paduapoint`](@ref) at position `(i, j)` is a Padua point.

# Examples
```jldoctest
julia> pointornothing(i, j, n) = PaduaTransforms.ispadua(i, j) ? PaduaTransforms.paduapoint(Float64, j, i, n) : nothing
pointornothing (generic function with 1 method)

julia> [pointornothing(y, x, 2) for y in 0:3, x in 0:2]
4×3 Matrix{Union{Nothing, Tuple{Float64, Float64}}}:
 (1.0, 1.0)   nothing      (-1.0, 1.0)
 nothing      (0.0, 0.5)   nothing
 (1.0, -0.5)  nothing      (-1.0, -0.5)
 nothing      (0.0, -1.0)  nothing
```

"""
ispadua(i, j) = iseven(i - j)



"""
    getpaduapoints([f::Function,], [T=Float64,] n)

evaluates the function `f` on the Padua points (of type `T`) for degree `n`.

If no function `f` is provided, `getpaduapoints` returns a matrix containing
the x and y components of the paduapoints as columns. If `f` returns a single value,
`getpaduapoints` returns a `Vector{T}`. If `f` returns a tuple or other iterable,
`getpaduapoints` returns a `Matrix{T}` where the i-th column represents the i-th
entry in `f` applied to all the Padua points.

# Examples
```jldoctest
julia> getpaduapoints(2)
6×2 Matrix{Float64}:
  1.0   1.0
  1.0  -0.5
  0.0   0.5
  0.0  -1.0
 -1.0   1.0
 -1.0  -0.5

julia> getpaduapoints(Float32, 1)
3×2 Matrix{Float32}:
  1.0   1.0
  1.0  -1.0
 -1.0   0.0

julia> getpaduapoints(Float32, 2) do x, y; x*y; end
6-element Vector{Float32}:
  1.0
 -0.50000006
  0.0
 -0.0
 -1.0
  0.50000006

julia> getpaduapoints(2) do x, y; x*y, 5*x*y; end
6×2 Matrix{Float64}:
  1.0   5.0
 -0.5  -2.5
  0.0   0.0
 -0.0  -0.0
 -1.0  -5.0
  0.5   2.5

```
"""
function getpaduapoints(f::Function, ::Type{T}, n) where T
    D = length(f(zero(T), zero(T)))
    out = Matrix{T}(undef, getpaduanum(n), D)

    i = 1
    for x in 0:n
        for y in 0:n+1
            if ispadua(x, y)
                px, py = paduapoint(T, x, y, n)
                out[i, :] .= f(px, py)
                i += 1
            end
        end
    end

    # Should return vector instead of matrix in 1D case
    return D == 1 ? vec(out) : out
end

getpaduapoints(f::Function, n) = getpaduapoints(f, Float64, n)
getpaduapoints(T, n) = getpaduapoints((x, y) -> (x, y), T, n)
getpaduapoints(n) = getpaduapoints(Float64, n)

## Padua Transform ##

struct PaduaTransformPlan{T, P}
    degree::Int
    vals::Matrix{T}
    dctplan::P
end

"""
    PaduaTransformPlan{T}(n::Integer)

create plan to compute coefficients of Chebyshev polynomials in 2D up to total degree `n`
using the Padua transform.

"""
function PaduaTransformPlan{T}(degree::Integer) where T
    vals = Matrix{T}(undef, degree + 2, degree + 1)
    plan = FFTW.plan_r2r!(vals, FFTW.REDFT00)

    PaduaTransformPlan{T, typeof(plan)}(degree, vals, plan)
end

"""
    weight!(mat::AbstractMatrix, degree::Integer)

weight fourier coefficients to obtain Chebyshev coefficients as part of a [`paduatransform!`](@ref).
The weighting factor applied to the coefficients is

```math
w = \\frac{1}{n(n+1)} ⋅ \\begin{cases}
    \\frac{1}{2} & \\textrm{if on vertex}   \\\\
    1           & \\textrm{if on edge}    \\\\
    2           & \\textrm{if in interior} \\\\
\\end{cases}
```

# Examples
```julia-repl
julia> weight!(ones(4+2, 4+1), 4)
6×5 Matrix{Float64}:
 0.025  0.05  0.05  0.05  0.025
 0.05   0.1   0.1   0.1   0.05
 0.05   0.1   0.1   0.1   0.05
 0.05   0.1   0.1   0.1   0.05
 0.05   0.1   0.1   0.1   0.05
 0.025  0.05  0.05  0.05  0.025
```
"""
function weight!(mat::AbstractMatrix{T}, degree::Integer) where T
    rmul!(mat, T(2) / T(degree * (degree + 1)))
    rmul!(@view(mat[1, :]), T(0.5))
    rmul!(@view(mat[end, :]), T(0.5))
    rmul!(@view(mat[:, 1]), T(0.5))
    rmul!(@view(mat[:, end]), T(0.5))

    mat
end

"""
    tovalsmat!(mat::Matrix, from::AbstractVector, degree::Integer)

write values of function evaluated at Padua points from `from` to matrix `mat`.

# Examples
```jldoctest
julia> PaduaTransforms.tovalsmat!(ones(3 + 2, 3 + 1), 1:getpaduanum(3), 3)
5×4 Matrix{Float64}:
 1.0  0.0  6.0   0.0
 0.0  4.0  0.0   9.0
 2.0  0.0  7.0   0.0
 0.0  5.0  0.0  10.0
 3.0  0.0  8.0   0.0

julia> PaduaTransforms.tovalsmat!(ones(2 + 2, 2 + 1), 1:getpaduanum(2), 2)
4×3 Matrix{Float64}:
 1.0  0.0  5.0
 0.0  3.0  0.0
 2.0  0.0  6.0
 0.0  4.0  0.0
```
"""
function tovalsmat!(mat::Matrix{T}, from::AbstractVector, degree::Integer) where T
    axes(from, 1) == 1:getpaduanum(degree) || throw(ArgumentError(
        "axes of from must be (1:getpaduanum(degree),)"))
    size(mat) == (degree + 2, degree + 1) || throw(ArgumentError(
        "axes of vals mat must be (degree + 2, degree + 1)"))

    if isodd(degree)
        # x 0
        # 0 x
        # x 0

        @inbounds for i in 1:length(from)
            mat[2i - 1] = from[i]
            mat[2i] = zero(T)
        end
    else
        @assert iseven(degree)
        # x 0 x
        # 0 x 0
        # x 0 x
        # 0 x 0

        valspercol = (degree + 2) ÷ 2

        # odd columns (j is column index)
        for j in 1:2:degree + 1, i in 1:valspercol
            k = (j - 1) * valspercol + i
            @inbounds mat[2i - 1, j] = from[k]
            @inbounds mat[2i, j] = zero(T)
        end

        # even columns
        for j in 2:2:degree + 1, i in 1:valspercol
            k = (j - 1) * valspercol + i
            @inbounds mat[2i - 1, j] = zero(T)
            @inbounds mat[2i, j] = from[k]
        end
    end

    mat
end

"""
    fromcoeffsmat!(to::AbstractVector, mat::Matrix, degree::Integer, ::Val{lex})

write Chebyshev coefficients from `mat` into vector `to`. `lex::Bool` determines whether
coefficients should be written in lexigraphical order or not. The lower right triangle does
not get written into `to`. These would represent greater polynomial degrees than `degree`.

If `lex` is `Val(true)` the coefficients correspond to the following basis polynomials

```math
T_0(x) T_0(y), T_1(x) T_0(y), T_0(x) T_1(y), T_2(x) T_0(y), T_1(x) T_1(y), T_0(x) T_2(y), ...
```

else if `lex` is `Val(false)` they correspond to

```math
T_0(x) T_0(y), T_0(x) T_1(y), T_1(x) T_0(y), T_0(x) T_2(y), T_1(x) T_1(y), T_2(x) T_0(y), ...
```

# Examples
```julia-repl
julia> mat = [(x, y) for y in 0:2+1, x in 0:2]
4×3 Matrix{Tuple{Int64, Int64}}:
 (0, 0)  (1, 0)  (2, 0)
 (0, 1)  (1, 1)  (2, 1)
 (0, 2)  (1, 2)  (2, 2)
 (0, 3)  (1, 3)  (2, 3)

julia> to1 = similar(mat, getpaduanum(2)); to2 = similar(mat, getpaduanum(2));

julia> fromcoeffsmat!(to1, mat, 2, Val(true))
6-element Vector{Tuple{Int64, Int64}}:
 (0, 0)
 (1, 0)
 (0, 1)
 (2, 0)
 (1, 1)
 (0, 2)

julia> fromcoeffsmat!(to2, mat, 2, Val(false))
6-element Vector{Tuple{Int64, Int64}}:
 (0, 0)
 (0, 1)
 (1, 0)
 (0, 2)
 (1, 1)
 (2, 0)
```
"""
function fromcoeffsmat!(to::AbstractVector, mat::Matrix, degree::Integer, ::Val{false})
    length(to) == getpaduanum(degree) || throw(ArgumentError(
        "length of output vector must equal paduanum(degree)"))
    axes(mat) == (1:(degree + 2), 1:(degree + 1)) || throw(ArgumentError(
        "axes of coeffs mat must be (1:(degree + 2), 1:(degree + 1))"))

    n = firstindex(to)
    for d in 1:degree + 1
        for ix in 1:d
            iy = d - ix + 1
            @assert ix + iy == d + 1 "ix and iy must lie on d-th diagonal"

            to[n] = mat[iy, ix]
            n += 1
        end
    end

    to
end

function fromcoeffsmat!(to::AbstractVector, mat::Matrix, degree::Integer, ::Val{true})
    length(to) == getpaduanum(degree) || throw(ArgumentError(
        "length of output must equal paduanum(degree)"))
    size(mat) == (degree + 2, degree + 1) || throw(ArgumentError(
        "axes of coeffs mat must be (1:(degree + 2), 1:(degree + 1))"))

    n = firstindex(to)
    for d in 1:degree + 1
        for iy in 1:d
            ix = d - iy + 1
            @assert ix + iy == d + 1 "ix and iy must lie on d-th diagonal"

            to[n] = mat[iy, ix]
            n += 1
        end
    end

    to
end

"""
    fromcoeffsmat!(to::AbstractMatrix, mat::Matrix, degree::Integer)

copy Chebyshev coefficients from `mat` to `to` without copying coefficients coresponding to
total degree greater than `degree`.

# Examples
```jldoctest
julia> PaduaTransforms.fromcoeffsmat!(zeros(4, 4), reshape(1:20, 5, 4), 3)
4×4 Matrix{Float64}:
 1.0  6.0  11.0  16.0
 2.0  7.0  12.0   0.0
 3.0  8.0   0.0   0.0
 4.0  0.0   0.0   0.0
```
"""
function fromcoeffsmat!(to::AbstractMatrix, mat::AbstractMatrix, degree::Integer)
    axes(to) == (1:degree + 1, 1:degree + 1) || throw(ArgumentError(
        "axes of output must be (1:(degree + 2), 1:(degree + 1))"))
    axes(mat) == (1:degree + 2, 1:degree + 1) || throw(ArgumentError(
        "axes of coeffs mat must be (1:(degree + 2), 1:(degree + 1))"))

    for j in 1:degree + 1
        for i in 1:degree + 2 - j
            @inbounds to[i, j] = mat[i, j]
        end
    end

    to
end

"""
    paduatransform!(coeffs, P::PaduaTransformPlan, vals[, lex])

obtain coefficients of Chebyshev polynomials on 2D via the Padua transform, given values
`vals` evaluated at the Padua points. Coefficients will be written into `coeffs`, which
should either be a matrix, a vector or an iterable of either.

If a matrix of `vals` and a corresponding iterable of `coeffs` matrixes or
vectors is given, each `vals` column will be transformed seperately in a multidimensional
Padua transform.

`lex` determines the order in which coefficients are written into `out` if `out` is a
vector. See [`fromcoeffsmat!`](@ref) for details.

!!! warning
    if `coeffs` is a matrix, make sure that all entries in the lower right diagonal are zero
    as these will not get overwritten.

# Examples
```jldoctest
julia> plan = PaduaTransformPlan{Float64}(2);

julia> f(x, y) = 3 + 4x + 5 * x*y, 6 + 7y
f (generic function with 1 method)

julia> vals = getpaduapoints(f, 2)
6×2 Matrix{Float64}:
 12.0  13.0
  4.5   2.5
  3.0   9.5
  3.0  -1.0
 -6.0  13.0
  1.5   2.5

julia> paduatransform!(zeros(3, 3), plan, vals[:, 1])
3×3 Matrix{Float64}:
 3.0  4.0  0.0
 0.0  5.0  0.0
 0.0  0.0  0.0

julia> paduatransform!(zeros(6), plan, vals[:, 2], Val(true))
6-element Vector{Float64}:
 6.0
 0.0
 7.0
 0.0
 0.0
 0.0

julia> paduatransform!((zeros(3, 3), zeros(3, 3)), plan, vals)
([3.0 4.0 0.0; 0.0 5.0 0.0; 0.0 0.0 0.0], [6.0 0.0 0.0; 7.0 0.0 0.0; 0.0 0.0 0.0])
```
"""
function paduatransform!(P::PaduaTransformPlan)
    coeffs = P.dctplan * P.vals
    weight!(coeffs, P.degree)

    coeffs
end

function paduatransform!(out, P::PaduaTransformPlan, vals::AbstractVector{<: Number}, args...)
    tovalsmat!(P.vals, vals, P.degree)
    coeffs = paduatransform!(P)
    fromcoeffsmat!(out, coeffs, P.degree, args...)
end

function paduatransform!(coeffs, P::PaduaTransformPlan, vals, args...)
    length(coeffs) == length(vals)|| throw(ArgumentError(
        "the length of coeffs and vals iterables must match"))

    # can't use eltype to check eltype beforehand because eachcol/eachrow are type unstable
    # and I want eachcol/eachrow to be usable as an iterable of vectors
    for (c, v::AbstractVector{<:Number}) in zip(coeffs, vals)
        paduatransform!(c, P, v, args...)
    end

    coeffs
end

function paduatransform!(coeffs, P::PaduaTransformPlan, vals::AbstractMatrix{<:Number}, args...)
    paduatransform!(coeffs, P, eachcol(vals), args...)
end

## Inverse Padua Transform ##

struct InvPaduaTransformPlan{T, P}
    degree::Int
    coeffs::Matrix{T}
    dctplan::P
end

"""
    InvPaduaTransformPlan{T}(n::Integer)

create plan to compute values on Padua points, given coefficients of Chebyshev polynomials
up to total degree `n`.
"""
function InvPaduaTransformPlan{T}(degree::Integer) where T
    coeffs = Matrix{T}(undef, degree + 2, degree + 1)
    iplan = FFTW.plan_r2r!(coeffs, FFTW.REDFT00)

    InvPaduaTransformPlan{T, typeof(iplan)}(degree, coeffs, iplan)
end

"""
    tocoeffsmat!(mat::AbstractMatrix, coeffs::AbstractMatrix)

writes coefficients in `coeffs` into matrix `mat` for the [`invpaduatransform!`](@ref).

# Examples
```jldoctest
julia> PaduaTransforms.tocoeffsmat!(zeros(5, 4), reshape(1:16, 4, 4))
5×4 Matrix{Float64}:
 1.0  5.0   9.0  13.0
 2.0  6.0  10.0  14.0
 3.0  7.0  11.0  15.0
 4.0  8.0  12.0  16.0
 0.0  0.0   0.0   0.0
```
"""
function tocoeffsmat!(mat::AbstractMatrix{T}, coeffs::AbstractMatrix) where T
    mat[1:end-1, :]  = coeffs
    mat[    end, :] .= zero(T)

    mat
end

"""
    invweight!(coeffs::AbstractMatrix)

weight Chebyshev coefficients before the Fourier transform for the [`invpaduatransform!`](@ref).
using the weighting

```math
w = \\begin{cases}
    1            & \\textrm{if on vertex}   \\\\
    \\frac{1}{2} & \\textrm{if on edge}     \\\\
    \\frac{1}{4} & \\textrm{if in interior} \\\\
\\end{cases}
```

# Examples
```jldoctest
julia> PaduaTransforms.invweight!(ones(5, 5))
5×5 Matrix{Float64}:
 1.0  0.5   0.5   0.5   1.0
 0.5  0.25  0.25  0.25  0.5
 0.5  0.25  0.25  0.25  0.5
 0.5  0.25  0.25  0.25  0.5
 1.0  0.5   0.5   0.5   1.0
```
"""
function invweight!(coeffs::AbstractMatrix{T}) where T
    rmul!(@view(coeffs[:,2:end-1]), T(0.5))
    rmul!(@view(coeffs[2:end-1, :]), T(0.5))

    coeffs
end

"""
    fromvalsmat!(to::AbstractVector, mat::AbstractMatrix, n::Integer)

write values from `mat` into the vector `to` after an [`invpaduatransform!`](@ref) of total
degree `n`.
"""
function fromvalsmat!(to::AbstractVector, mat::AbstractMatrix, degree::Integer)
    axes(to, 1) == 1:getpaduanum(degree) || throw(Argumenterror(
        "axes of output must be (1:getpaduanum(degree),)"))
    axes(mat) == (1:degree + 2, 1:degree + 1) || throw(Argumenterror(
        "axes of vals mat must be (1:degree + 2, 1:degree + 1)"))

    if isodd(degree)
        # x 0
        # 0 x
        # x 0

        @inbounds for i in 1:length(to)
            to[i] = mat[2i - 1]
        end
    else
        @assert iseven(degree)
        # x 0 x
        # 0 x 0
        # x 0 x
        # 0 x 0

        valspercol = (degree + 2) ÷ 2

        # odd columns (j is column index)
        for j in 1:2:degree + 1, i in 1:valspercol
            k = (j - 1) * valspercol + i
            @inbounds to[k] = mat[2i - 1, j]
        end

        # even columns
        for j in 2:2:degree + 1, i in 1:valspercol
            k = (j - 1) * valspercol + i
            @inbounds to[k] = mat[2i, j]
        end
    end

    to
end

function invpaduatransform!(IP::InvPaduaTransformPlan)
    invweight!(IP.coeffs)
    IP.dctplan * IP.coeffs

    IP.coeffs
end

"""
    invpaduatransform!(vals::AbstractVector, IP::InvPaduaTransformPlan, coeffs::AbstractMatrix)

evaluates the polynomial defined by the coefficients of Chebyshev polynomials `coeffs` on
the Padua points using the inverse transform plan `IP` and writes the resulting values into
`vals`.

If a `vals` matrix and a corresponding iterable of `coeffs` matrixes is given,
each `coeffs` matrix will be transformed seperately into a column in `vals`.

# Examples
```jldoctest
julia> iplan = InvPaduaTransformPlan{Float64}(2);

julia> coeffs = [3 4 0; 0 5 0; 0 0 0]
3×3 Matrix{Int64}:
 3  4  0
 0  5  0
 0  0  0

julia> invpaduatransform!(zeros(6), iplan, coeffs)
6-element Vector{Float64}:
 12.0
  4.5
  3.0
  3.0
 -6.0
  1.5

julia> getpaduapoints(2) do x, y; 3 + 4x + 5 * x*y; end
6-element Vector{Float64}:
 12.0
  4.5
  3.0
  3.0
 -6.0
  1.4999999999999996
```
"""
function invpaduatransform!(vals::AbstractVector{<:Number}, IP::InvPaduaTransformPlan, coeffs)
    tocoeffsmat!(IP.coeffs, coeffs)
    invpaduatransform!(IP)
    fromvalsmat!(vals, IP.coeffs, IP.degree)
end

function invpaduatransform!(vals, IP::InvPaduaTransformPlan, coeffs)
    length(coeffs) == length(vals)|| throw(ArgumentError(
        "the length of coeffs and vals iterables must match"))

    for (v::AbstractVector{<:Number}, c) in zip(vals, coeffs)
        invpaduatransform!(v, IP, c)
    end

    vals
end

function invpaduatransform!(vals::AbstractMatrix{<:Number}, IP::InvPaduaTransformPlan, coeffs, args...)
    invpaduatransform!(eachcol(vals), IP, coeffs, args...)
end

end  # PaduaTransforms
