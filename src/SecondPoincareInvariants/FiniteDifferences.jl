module FiniteDifferences

import ...PoincareInvariants: compute!, getpoints, getpointnum, getpointspec
import ..SecondPoincareInvariants: SecondPoincareInvariant

using ...PoincareInvariants: @argcheck

using LinearAlgebra: dot
using Base: Callable
using StaticArrays: SVector

struct FiniteDiffPlan{T, D}
    ∂x::NTuple{D, Vector{T}}
    ∂y::NTuple{D, Vector{T}}
    simpweights::Vector{T}
end

function FiniteDiffPlan{T, D}(Ω, ps::NTuple{2, Int}) where {T, D}
    nx, ny = ps
    N = nx * ny

    ∂x = ntuple(_ -> Vector{T}(undef, N), D)
    ∂y = ntuple(_ -> Vector{T}(undef, N), D)
    simpweights = getsimpweights(T, nx, ny)

    FiniteDiffPlan{T, D}(∂x, ∂y, simpweights)
end

function compute!(
    pinv::SecondPoincareInvariant{T, D, ΩT, <:Any, P}, phasepoints, t, p
) where {T, D, ΩT<:Callable, P <: FiniteDiffPlan}
    nx, ny = pinv.pointspec
    N = nx * ny

    differentiate!(pinv.plan.∂x, pinv.plan.∂y, phasepoints, (nx, ny))

    I = zero(T)
    w = pinv.plan.simpweights

    @inbounds for i in 1:N
        pnti = SVector{D}(ntuple(d -> phasepoints[d][i], D))
        ∂xi = SVector{D}(ntuple(d -> pinv.plan.∂x[d][i], D))
        ∂yi = SVector{D}(ntuple(d -> pinv.plan.∂y[d][i], D))

        I += w[i] * dot(∂yi, pinv.Ω(pnti, t, p), ∂xi)
    end

    return I
end

## getpoints, getpointnum and getpointspec ##

getpointnum((nx, ny)::NTuple{2, Integer}) = nx * ny

nextodd(n::Integer)::Int = n <= 3 ? 3 : 2 * fld(n, 2) + 1

function getpointspec(dims::NTuple{2, Integer}, ::Type{<:FiniteDiffPlan})
    return (nextodd(dims[1]), nextodd(dims[2]))
end

function getpointspec(N::Integer, ::Type{<:FiniteDiffPlan})
    n = nextodd(ceil(Int, sqrt(Int(N))))
    return (n, n)
end

function getpoints(f, ::Type{T}, dims::NTuple{2, Integer}, ::Type{<:FiniteDiffPlan}) where T
    D = length(f(zero(T), zero(T)))
    nx, ny = dims
    N = nx * ny
    out = ntuple(_ -> Vector{T}(undef, N), D)

    i = 1
    for x in range(0, 1, length=nx)
        for y in range(0, 1, length=ny)
            fpnt = f(x, y)
            for d in 1:D
                out[d][i] = fpnt[d]
            end
            i += 1
        end
    end

    # Should return vector [...] instead of ([...],) in 1D case
    return D == 1 ? out[1] : out
end

## Differentiation ##

# convert a cartesian index (x, y) to a linear index and index into vector
_tolin(x::Integer, y::Integer, ny::Integer) = (x - 1) * ny + y

# evaluate derivative of quadratic approxmiation around (x0, y0) at (x0 + x, y0 + y)
# We have approximation f(x, y) = c + cx*x + cy*y + cxx*x^2 + cxy*x*y + cyy*y^2
function _diff(
    vals::AbstractVector{T}, x0::Integer, y0::Integer,
    x::Integer, y::Integer, nx::Integer, ny::Integer,
    ::Val{edge}
) where {T, edge}
    f(ix, iy) = @inbounds vals[_tolin(x0 + ix, y0 + iy, ny)]

    # All coefficients are multiplied by 2*Δx or 2*Δy
    bx = f(1, 0); ax = f(-1, 0)
    fx = cx = bx - ax

    by = f(0, 1); ay = f(0, -1)
    fy = cy = by - ay

    if edge  # meaning if (x, y) not (0, 0)
        f00 = f(0, 0)
        cxx = bx + ax - 2*f00
        cyy = by + ay - 2*f00
        cxy = (f(1, 1) + f(-1, -1) - f(1, -1) - f(-1, 1)) / 2

        fx += 2*cxx * x + cxy * y
        fy += 2*cyy * y + cxy * x
    end

    invΔx = nx - 1; invΔy = ny - 1
    return fx * invΔx / 2, fy * invΔy / 2
end

function differentiate!(
    ∂x::AbstractVector, ∂y::AbstractVector, vals::AbstractVector{T},
    dims::NTuple{2, Integer}
) where T
    nx, ny = dims
    Base.require_one_based_indexing(∂x, ∂y, vals)
    @argcheck nx >= 3 && ny >= 3 "must have at least 3×3 points to calculate both derivatives"
    @argcheck axes(∂x, 1) == axes(∂y, 1) == axes(vals, 1) "axes of derivatives and " *
        "values must match"
    @argcheck length(vals) == nx * ny "length of vals and length implied by dims must match"

    _li(x, y) = _tolin(x, y, ny)

    # corners
    ∂x[_li(1 ,  1)], ∂y[_li(1 ,  1)] = _diff(vals,    2,    2, -1, -1, nx, ny, Val(true))
    ∂x[_li(nx,  1)], ∂y[_li(nx,  1)] = _diff(vals, nx-1,    2, +1, -1, nx, ny, Val(true))
    ∂x[_li(nx, ny)], ∂y[_li(nx, ny)] = _diff(vals, nx-1, ny-1, +1, +1, nx, ny, Val(true))
    ∂x[_li(1 , ny)], ∂y[_li(1 , ny)] = _diff(vals,    2, ny-1, -1, +1, nx, ny, Val(true))

    # edges
    @inbounds for y in 2:ny-1
        ∂x[_li(1 ,  y)], ∂y[_li(1 ,  y)] = _diff(vals,    2, y, -1, 0, nx, ny, Val(true))
        ∂x[_li(nx,  y)], ∂y[_li(nx,  y)] = _diff(vals, nx-1, y, +1, 0, nx, ny, Val(true))
    end

    @inbounds for x in 2:nx-1
        ∂x[_li(x,  1)], ∂y[_li(x,  1)] = _diff(vals, x,    2, 0, -1, nx, ny, Val(true))
        ∂x[_li(x, ny)], ∂y[_li(x, ny)] = _diff(vals, x, ny-1, 0, +1, nx, ny, Val(true))
    end

    # interior (uses central differences)
    @inbounds for x in 2:nx-1
        for y in 2:ny-1
            ∂x[_li(x, y)], ∂y[_li(x, y)] = _diff(vals, x, y, 0, 0, nx, ny, Val(false))
        end
    end

    return ∂x, ∂y
end

function differentiate!(∂x, ∂y, vals, dims::NTuple{2, Integer})
    @argcheck length(∂x) == length(∂y) == length(vals) "number of derivative and vals " *
        "vectors must match"

    for (∂xi, ∂yi, valsi) in zip(∂x, ∂y, vals)
        differentiate!(∂xi, ∂yi, valsi, dims)
    end

    return ∂x, ∂y
end

## Integration ##

function _sw(i, n)
    if i == 1 || i == n
        return 1
    elseif iseven(i)
        return 4
    else  # odd
        return 2
    end
end

function getsimpweights(::Type{T}, nx, ny) where T
    w = Matrix{T}(undef, ny, nx)
    for x in 1:nx
        wx = _sw(x, nx)
        for y in 1:ny
            wy = _sw(y, ny)
            w[y, x] = T(wx * wy) / (9 * (nx - 1) * (ny - 1))
        end
    end
    return vec(w)
end

end  # FiniteDifferences
