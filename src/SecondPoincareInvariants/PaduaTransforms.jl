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
up to total degree `n`. This number is equal to the number of coefficients.

# Examples
```julia-repl
julia> getpaduanum(13)
105
```
"""
getpaduanum(degree) = (degree + 1) * (degree + 2) ÷ 2

"""
    getdegree(N)

calculates total degree of corresponding polynomial, given the number of coefficients or
Padua points `N`.

# Examples
```julia-repl
julia> getdegree(105)
13
```
"""
function getdegree(pointnum)
	d = (sqrt(1 + 8pointnum) - 3) / 2
	isinteger(d) ? Int(d) : throw(ArgumentError(
		"number of Padua points or coeffs must be (n + 1) * (n + 2) ÷ 2"))
end

"""
    nextpaduanum(N)

get next valid paduanum ≥ `N`.

# Examples

```julia-repl
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
    chebyshevpoint(T, ix, iy, nx, ny)

return Chebyshev point at grid position `(ix, iy)` on an `nx × ny` grid as an `SVector{2, T}`.

`z = (cos(π * ix / nx), cos(π * iy / ny))`

# Examples
```julia-repl
julia> [chebyshevpoint(Float64, x, y, 2, 3) for y in 0:3, x in 0:2]
4×3 Matrix{SVector{2, Float64}}:
 [1.0, 1.0]   [0.0, 1.0]   [-1.0, 1.0]
 [1.0, 0.5]   [0.0, 0.5]   [-1.0, 0.5]
 [1.0, -0.5]  [0.0, -0.5]  [-1.0, -0.5]
 [1.0, -1.0]  [0.0, -1.0]  [-1.0, -1.0]
```
"""
function chebyshevpoint(::Type{T}, ix, iy, nx, ny) where T
    x = cospi(T(ix) / T(nx))
    y = cospi(T(iy) / T(ny))
    SVector{2, T}(x, y)
end

"""
    paduapoint(T, ix, iy, n)

returns [`chebyshevpoint`](@ref) at position `(ix, iy)` on an `n × (n + 1)` grid.

Only every second point is actually a Padua point, though. Check with [`ispadua`](@ref).

`z = (cos(π * ix / n), cos(π * iy / (n + 1)))`

# Examples
```julia-repl
julia> [paduapoint(Float32, x, y, 1) for y in 0:1+1, x in 0:1]
3×2 Matrix{SVector{2, Float32}}:
 [1.0, 1.0]   [-1.0, 1.0]
 [1.0, 0.0]   [-1.0, 0.0]
 [1.0, -1.0]  [-1.0, -1.0]
```

"""
function paduapoint(::Type{T}, ix, iy, n) where T
    chebyshevpoint(T, ix, iy, n, n + 1)
end

# This checks whether a point on the chebyshev grid at position (i, j) is also a Padua point
"""
    ispadua(ix, iy)

checks whether [`chebyshevpoint`](@ref) at position `(ix, iy)` is a Padua point.

# Examples
```julia-repl
julia> [ispadua(x, y) ? chebyshevpoint(Float64, x, y, 2, 3) : nothing for y in 0:3, x in 0:2]
4×3 Matrix{Union{Nothing, SVector{2, Float64}}}:
 [1.0, 1.0]   nothing      [-1.0, 1.0]
 nothing      [0.0, 0.5]   nothing
 [1.0, -0.5]  nothing      [-1.0, -0.5]
 nothing      [0.0, -1.0]  nothing
```

"""
ispadua(i, j) = iseven(i - j)

"""
    getpaduapoints([T=Float64, ]n)

returns `(n + 1)(n + 2) / 2` Padua points.

`getpaduapoints(n)` is equivalent to

```
[paduapoint(Float64, x, y, n) for y in 0:n+1, x in 0:n if ispadua(x, y)]
```

# Examples

```julia-repl
julia> getpaduapoints(Float64, 1)
3-element Vector{SVector{2, Float64}}:
 [1.0, 1.0]
 [1.0, -1.0]
 [-1.0, 0.0]
```

"""
getpaduapoints(::Type{T}, n) where T = [paduapoint(T, x, y, n) for y in 0:n+1, x in 0:n if ispadua(x, y)]
getpaduapoints(n) = getpaduapoints(Float64, n)

## Padua Transform ##

struct PaduaTransformPlan{T, P}
    degree::Int
    vals::Matrix{T}
    dctplan::P
end

"""
    PaduaTransformPlan{T}(degree::Integer)

create plan to compute coefficients of Chebyshev polynomials in 2D up to degree `degree`
using the Padua transform, as implemented in [this](https://link.springer.com/article/10.1007/s11075-010-9373-1)
paper.

See also: [`paduatransform!`](@ref)
"""
function PaduaTransformPlan{T}(degree::Integer) where T
    vals = Matrix{T}(undef, degree + 2, degree + 1)
    plan = FFTW.plan_r2r!(vals, FFTW.REDFT00)

    PaduaTransformPlan{T, typeof(plan)}(degree, vals, plan)
end

"""
    weight!(mat::AbstractMatrix, degree::Integer)

weight fourier coefficients to obtain Chebyshev coefficients.

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
    tovalsmat!(mat::Matrix, vec::AbstractVector, degree::Integer)

write values of function evaluated at Padua points from `vec` to to matrix `mat`.

# Examples
```julia-repl
julia> tovalsmat!(ones(3 + 2, 3 + 1), 1:getpaduanum(3), 3)
5×4 Matrix{Float64}:
 1.0  0.0  6.0   0.0
 0.0  4.0  0.0   9.0
 2.0  0.0  7.0   0.0
 0.0  5.0  0.0  10.0
 3.0  0.0  8.0   0.0

julia> tovalsmat!(ones(2 + 2, 2 + 1), 1:getpaduanum(2), 2)
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

        mat[1:2:end] .= from
		mat[2:2:end] .= zero(T)
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
    fromcoeffsmat!(to::AbstractVector, mat::AbstractMatrix, degree::Integer, ::Val{lex})

write Chebyshev coefficients from `mat` into vector `to`. `lex::Bool` determines whether
coefficients should be written in lexigraphical order or not. (See examples)

The lower right triangle does not get written into `to`. These would represent higher
polynomial degrees than `degree`.

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
    fromcoeffsmat!(to::AbstractMatrix, mat::AbstractMatrix, degree::Integer)

copy Chebyshev coefficients from `mat` to `to` without copying coefficients coresponding to
total degree > `degree`.

# Examples
```
julia> fromcoeffsmat!(zeros(5, 4), reshape(1:20, 5, 4), 3)
5×4 Matrix{Float64}:
 1.0  6.0  11.0  16.0
 2.0  7.0  12.0   0.0
 3.0  8.0   0.0   0.0
 4.0  0.0   0.0   0.0
 0.0  0.0   0.0   0.0
```
"""
function fromcoeffsmat!(to::AbstractMatrix, mat::AbstractMatrix, degree::Integer)
	axes(to) == axes(mat) || error()
	size(mat) == (degree + 2, degree + 1) || error()

	for j in 1:degree + 1
		for i in 1:degree + 2 - j
			@inbounds to[i, j] = mat[i, j]
		end
	end

	to
end

"""
    paduatransform!([out, ]P::PaduaTransformPlan, vals, args...)

obtain coefficients of Chebyshev polynomials on 2D via the Padua transform, given values
`vals` evaluated at the Padua points. Coefficients will either be returned as a matrix, or,
if `out` is given, written into `out`. `args` are passed to [`fromcoeffsmat!`](@ref).

# Examples
```julia-repl
julia> plan = PaduaTransformPlan{Float64}(3); points = getpaduapoints(3);

julia> g(v) = 3 + 4v[1] + 5 * v[2] * (2v[1]^2 - 1); vals = g.(points);

julia> paduatransform!(plan, vals)
5×4 Matrix{Float64}:
 3.0  4.0  0.0  0.0
 0.0  0.0  5.0  0.0
 0.0  0.0  0.0  0.0
 0.0  5.0  0.0  0.0
 0.0  0.0  4.0  3.0

julia> paduatransform!(zeros(getpaduanum(3)), plan, vals, Val(true))
[3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0]

julia> paduatransform!(zeros(5, 4), plan, vals)
5×4 Matrix{Float64}:
 3.0  4.0  0.0  0.0
 0.0  0.0  5.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
```
"""
function paduatransform!(P::PaduaTransformPlan)
    coeffs = P.dctplan * P.vals
    weight!(coeffs, P.degree)

    coeffs
end

function paduatransform!(P::PaduaTransformPlan, vals)
    tovalsmat!(P.vals, vals, P.degree)
    paduatransform!(P)
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
up to degree `n`.

See also: [`invpaduatransform!`](@ref)
"""
function InvPaduaTransformPlan{T}(degree::Integer) where T
    coeffs = Matrix{T}(undef, degree + 2, degree + 1)
    iplan = FFTW.plan_r2r!(coeffs, FFTW.REDFT00)

    InvPaduaTransformPlan{T, typeof(iplan)}(degree, coeffs, iplan)
end

"""
    invweight!(coeffs::AbstractMatrix)

weight Chebyshev coefficients before the Fourier transform in the inverse Padua transform.

# Examples
```julia-repl
julia> invweight!(ones(5, 5))
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


function fromvalsmat!(to::AbstractVector, mat::Matrix, degree::Integer)
	axes(to, 1) == 1:getpaduanum(degree) || error()
	size(mat) == (degree + 2, degree + 1) || error()

	if isodd(degree)
		# x 0
        # 0 x
        # x 0

		to .= mat[1:2:end]
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

function invpaduatransform!(IP::InvPaduaTransformPlan, coeffs::AbstractMatrix)
	IP.coeffs[:] = coeffs
	invpaduatransform!(IP)
end

function invpaduatransform!(vals::AbstractVector, IP::InvPaduaTransformPlan, coeffs::AbstractMatrix)
	invpaduatransform!(IP, coeffs)
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
