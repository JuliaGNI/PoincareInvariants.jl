"""
    PaduaTransforms

an implementation of the Padua transform and its inverse via the fast Fourier transform.
"""
module PaduaTransforms

using StaticArrays: SVector
using FFTW
using LinearAlgebra: rmul!

export getpaduanum, getdegree, nextpaduanum
export getpaduapoints
export PaduaTransformPlan, paduatransform!
export InvPaduaTransformPlan, invpaduatransform!

## Number of Padua Points and Degree ##

"""
    getpaduanum(n)

calculates number of Padua points needed to approximate a function using Chebyshev polynomials
up to total degree `n`. This number is equal to the number of coefficients. The formula is

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
    nextpaduanum(N)

get next valid number of Padua points ≥ `N`.

# Examples

```jldoctest
julia> nextpaduanum(104)
105
```
"""
function nextpaduanum(N)
    d = Int(cld(sqrt(1 + 8N) - 3, 2))
    getpaduanum(d)
end

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
    getpaduapoints([T=Float64,] n)

returns the Padua points

```math
\\textrm{Pad}_n = \\{(\\cos{\\frac{jπ}{n}}, \\cos{\\frac{iπ}{n + 1}}) \\; | \\;
    0 ≤ j ≤ n, \\; 0 ≤ i ≤ n + 1, \\; i - j \\; \\textrm{even} \\}
```

where each row is a point

# Examples

```jldoctest
julia> getpaduapoints(Float32, 1)
3×2 Matrix{Float32}:
  1.0   1.0
  1.0  -1.0
 -1.0   0.0
```
"""
function getpaduapoints(::Type{T}, n) where T
    out = Matrix{T}(undef, getpaduanum(n), 2)

    i = 1
    for x in 0:n
        for y in 0:n+1
            if ispadua(x, y)
                out[i, 1:2] .= paduapoint(T, x, y, n)
                i += 1
            end
        end
    end

    return out
end

getpaduapoints(n) = getpaduapoints(Float64, n)

"""
    getpaduapoints(f::Function, [T=Float64,] n)

evaluates the function `f` on the Padua points for degree `n`. If `f` returns a single value,
`getpaduapoints` returns a `Vector{T}`, else if `f` returns a tuple or other iterable
`getpaduapoints` returns a `Matrix{T}` where each row represents `f` applied to a Padua point.

# Examples
```jldoctest
julia> getpaduapoints(Float64, 2) do x, y; 3x - y, y^2; end
6×2 Matrix{Float64}:
  2.0  1.0
  3.5  0.25
 -0.5  0.25
  1.0  1.0
 -4.0  1.0
 -2.5  0.25
```
"""
function getpaduapoints(f::Function, ::Type{T}, n) where T
    D = length(f(zero(T), zero(T)))
    if D == 1
        out = Vector{T}(undef, getpaduanum(n))
    else
        out = Matrix{T}(undef, getpaduanum(n), D)
    end

    # Function barrier to alleviate issues from type instability
    _fillpoints!(out, f, T, n)

    return out
end

getpaduapoints(f::Function, n) = getpaduapoints(f, Float64, n)

function _fillpoints!(out::AbstractMatrix, f, ::Type{T}, n) where T
    i = 1
    for x in 0:n
        for y in 0:n+1
            if ispadua(x, y)
                v = paduapoint(T, x, y, n)
                out[i, :] .= f(v[1], v[2])
                i += 1
            end
        end
    end

    out
end

function _fillpoints!(out::AbstractVector, f, ::Type{T}, n) where T
    i = 1
    for x in 0:n
        for y in 0:n+1
            if ispadua(x, y)
                v = paduapoint(T, x, y, n)
                out[i] = f(v[1], v[2])
                i += 1
            end
        end
    end

    out
end

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
    rmul!(mat, T(2 / ( degree * (degree + 1) )))
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
    axes(from, 1) == 1:getpaduanum(degree) || error()
    size(mat) == (degree + 2, degree + 1) || error()

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
    length(to) == getpaduanum(degree) || error()
    axes(mat) == (1:(degree + 2), 1:(degree + 1)) || error()

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
    length(to) == getpaduanum(degree) || error()
    size(mat) == (degree + 2, degree + 1) || error()

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
    axes(to) == (1:degree + 1, 1:degree + 1) || error()
    axes(mat) == (1:degree + 2, 1:degree + 1) || error()

    for j in 1:degree + 1
        for i in 1:degree + 2 - j
            @inbounds to[i, j] = mat[i, j]
        end
    end

    to
end

"""
    paduatransform!(out, P::PaduaTransformPlan, vals[, lex])

obtain coefficients of Chebyshev polynomials on 2D via the Padua transform, given values
`vals` evaluated at the Padua points. Coefficients will be written into `out`, which should
either be a matrix or a vector.

if `out` is a matrix, make sure that all entries in the lower right diagonal are zero as
these will not get overwritten.

`lex` determines the order in which coefficeints are written into `out` if `out` is a vector.
See [`fromcoeffsmat!`](@ref) for details.

# Examples
```jldoctest
julia> plan = PaduaTransformPlan{Float64}(3);

julia> f(x, y) = 3 + 4x + 5 * y * (2x^2 - 1)
f (generic function with 1 method)

julia> vals = getpaduapoints(f, 3)
10-element Vector{Float64}:
 12.0
  7.0
  2.0
  3.232233047033631
  6.767766952966369
 -1.5000000000000004
  1.0000000000000004
  3.5000000000000013
  2.5355339059327378
 -4.535533905932738

julia> paduatransform!(zeros(4, 4), plan, vals)
4×4 Matrix{Float64}:
  3.0           4.0         0.0  1.4803e-16
 -5.92119e-16  -1.4803e-16  5.0  0.0
  0.0           0.0         0.0  0.0
 -2.96059e-16   0.0         0.0  0.0

julia> paduatransform!(zeros(getpaduanum(3)), plan, vals, Val(true))
10-element Vector{Float64}:
  3.0
  4.0
 -5.921189464667501e-16
  0.0
 -1.4802973661668753e-16
  0.0
  1.4802973661668753e-16
  5.0
  0.0
 -2.9605947323337506e-16

julia> paduatransform!(zeros(getpaduanum(3)), plan, vals, Val(false))
10-element Vector{Float64}:
  3.0
 -5.921189464667501e-16
  4.0
  0.0
 -1.4802973661668753e-16
  0.0
 -2.9605947323337506e-16
  0.0
  5.0
  1.4802973661668753e-16
```
"""
function paduatransform!(P::PaduaTransformPlan)
    coeffs = P.dctplan * P.vals
    weight!(coeffs, P.degree)

    coeffs
end

function paduatransform!(out, P::PaduaTransformPlan, vals, args...)
    tovalsmat!(P.vals, vals, P.degree)
    coeffs = paduatransform!(P)
    fromcoeffsmat!(out, coeffs, P.degree, args...)
end

"""
    paduatransform!(out::AbstractArray{<:Any, 3}, P::PaduaTransformPlan, vals::AbstractMatrix, args...)

transforms each column in `vals` and writes the resulting coefficients in a slice of `out`.
"""
function paduatransform!(out::AbstractArray{<:Any, 3}, P::PaduaTransformPlan, vals::AbstractMatrix, args...)
    axes(out, 3) == axes(vals, 2)|| error()

    @views for i in axes(out, 3)
        paduatransform!(out[:, :, i], P, vals[:, i], args...)
    end

    out
end

"""
    paduatransform!(out::Array{<:Any, 3}, P::PaduaTransformPlan, vals::AbstractVector{<:AbstractVector{T}}, args...)

transforms vector of vectors to Chebyshev coefficients and writes the resulting coefficients
in a slice of `out`. Each vector in the vector of vectors represents a point. Each slice of
`out` represents a transform of one dimension.
"""
function paduatransform!(out::Array{<:Any, 3}, P::PaduaTransformPlan, vals::AbstractVector{<:AbstractVector{T}}, args...) where T
    # Here, each column is a point and each row represents one dimension
    r = reinterpret(reshape, T, vals)
    axes(out, 3) == axes(r, 1) || error()

    @views for i in axes(out, 3)
        paduatransform!(out[:, :, i], P, r[i, :], args...)
    end

    out
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
    axes(to, 1) == 1:getpaduanum(degree) || error()
    axes(mat) == (1:degree + 2, 1:degree + 1) || error()

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

evaluates the polynomial defined by the coefficients of Chebyshev polynomials `coeffs` on the
Padua points using the inverse transform plan `IP` and writes the resulting values into `vals`.
"""
function invpaduatransform!(vals::AbstractVector, IP::InvPaduaTransformPlan, coeffs::AbstractMatrix)
    tocoeffsmat!(IP.coeffs, coeffs)
    invpaduatransform!(IP)
    fromvalsmat!(vals, IP.coeffs, IP.degree)
end

function invpaduatransform!(vals::AbstractMatrix, IP::InvPaduaTransformPlan, coeffs::AbstractArray{<:Any, 3})
    axes(coeffs, 3) == axes(vals, 2) || error()

    @views for i in axes(coeffs, 3)
        invpaduatransform!(vals[:, i], IP, coeffs[:, :, i])
    end

    vals
end

end  # Padua
