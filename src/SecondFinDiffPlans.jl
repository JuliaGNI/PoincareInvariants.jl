module SecondFinDiffPlans

import ..PoincareInvariants: compute!, getpoints, getpointnum, getpointspec
import ..PoincareInvariants: SecondPoincareInvariant

using ..PoincareInvariants: @argcheck

using LinearAlgebra: dot
using Base: Callable

struct SecondFinDiffPlan{T, D} end

function SecondFinDiffPlan{T, D}(ω, ps::NTuple{2, Int}) where {T, D}
    SecondFinDiffPlan{T, D}()
end

function compute!(
    pinv::SecondPoincareInvariant{T, D, ωT, <:Any, P}, vals, t, p
) where {T, D, ωT, P <: SecondFinDiffPlan}
    nx, ny = pinv.pointspec
    @argcheck size(vals) == (nx * ny, D) "Expected points mtarix to have size $((nx * ny, D))"

    I = zero(T)

    ∂xi = Vector{T}(undef, D)
    ∂yi = Vector{T}(undef, D)

    colviews = [view(vals, :, d) for d in 1:D]

    i = 1
    for ix in 1:nx
        for iy in 1:ny
            for d in 1:D
                ∂xi[d], ∂yi[d] = differentiate(colviews[d], ix, iy, (nx, ny))
            end

            # This if statement should hopefully get optimised away by the compiler
            if ωT <: Callable
                pnti = view(vals, i, :)
                ωi = pinv.ω(pnti, t, p)
            elseif ωT <: AbstractMatrix
                ωi = pinv.ω
            end

            w = getsimpweight(T, ix, iy, (nx, ny))

            I += w * dot(∂yi, ωi, ∂xi)

            i += 1
        end
    end

    return I
end

## getpoints, getpointnum and getpointspec ##

getpointnum((nx, ny)::NTuple{2, Integer}) = nx * ny

nextodd(n::Integer)::Int = n <= 3 ? 3 : 2 * fld(n, 2) + 1

function getpointspec(dims::NTuple{2, Integer}, ::Type{<:SecondFinDiffPlan})
    return (nextodd(dims[1]), nextodd(dims[2]))
end

function getpointspec(N::Integer, ::Type{<:SecondFinDiffPlan})
    n = nextodd(ceil(Int, sqrt(Int(N))))
    return (n, n)
end

function getpoints(f, ::Type{T}, dims::NTuple{2, Integer}, ::Type{<:SecondFinDiffPlan}) where T
    D = length(f(zero(T), zero(T)))
    nx, ny = dims
    N = nx * ny
    out = Matrix{T}(undef, N, D)

    i = 1
    for x in range(0, 1, length=nx)
        for y in range(0, 1, length=ny)
            out[i, :] .= f(x, y)
            i += 1
        end
    end

    # Should return vector instead of matrix in 1D case
    return D == 1 ? vec(out) : out
end

## Differentiation ##

function _getmid(ix, iy, nx, ny)
    x0, y0 = ix, iy
    edge = false

    if ix == 1
        x0 = 2
        edge = true
    end

    if iy == 1
        y0 = 2
        edge = true
    end

    if ix == nx
        x0 = nx - 1
        edge = true
    end

    if iy == ny
        y0 = ny - 1
        edge = true
    end

    return x0, y0, edge
end

# evaluate derivative of quadratic approxmiation around (x0, y0) at (x0 + x, y0 + y)
# We have approximation f(x, y) = c + cx*x + cy*y + cxx*x^2 + cxy*x*y + cyy*y^2
function differentiate(
    vals::AbstractVector{T}, ix::Integer, iy::Integer, (nx, ny)::NTuple{2, Integer}
) where T
    x0, y0, edge = _getmid(ix, iy, nx, ny)

    mat = reshape(vals, ny, nx)
    f(x, y) = @inbounds mat[y0 + y, x0 + x]

    # All coefficients are multiplied by 2*Δx or 2*Δy
    bx = f(1, 0); ax = f(-1, 0)
    fx = bx - ax  # cx

    by = f(0, 1); ay = f(0, -1)
    fy = by - ay  # cy

    if edge
        f00 = f(0, 0)
        cxx = bx + ax - 2*f00
        cyy = by + ay - 2*f00
        cxy = (f(1, 1) + f(-1, -1) - f(1, -1) - f(-1, 1)) / 2

        dx, dy = ix - x0, iy - y0
        fx += 2*cxx * dx + cxy * dy
        fy += 2*cyy * dy + cxy * dx
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
    else
        @assert isodd(i)
        return 2
    end
end

function getsimpweight(::Type{T}, x, y, (nx, ny)) where T
    return T(_sw(x, nx) * _sw(y, ny)) / (9 * (nx - 1) * (ny - 1))
end

end  # SecondFinDiffPlans
