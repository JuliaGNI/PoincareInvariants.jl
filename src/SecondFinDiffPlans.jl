module SecondFinDiffPlans

import ..PoincareInvariants: compute!, getpoints, getpointnum, getpointspec
import ..PoincareInvariants: SecondPoincareInvariant

using ..PoincareInvariants: @argcheck

using LinearAlgebra: dot

struct SecondFinDiffPlan{T, D} end

function SecondFinDiffPlan{T, D}(ω, ps::NTuple{2, Int}) where {T, D}
    SecondFinDiffPlan{T, D}()
end

function compute!(
    pinv::SecondPoincareInvariant{T, D, ωT, <:Any, P}, t::Real, p
) where {T, D, ωT, P <: SecondFinDiffPlan}
    nx, ny = pinv.pointspec
    points = pinv.points

    @argcheck size(points) == (nx * ny, D) "Expected points mtarix to have size $((nx * ny, D))"

    I = zero(T)

    ∂xi = Vector{T}(undef, D)
    ∂yi = Vector{T}(undef, D)

    colviews = [view(points, :, d) for d in 1:D]

    i = 1
    for ix in 1:nx
        for iy in 1:ny
            for d in 1:D
                ∂xi[d], ∂yi[d] = differentiate(colviews[d], ix, iy, (nx, ny))
            end

            # This if statement should hopefully get optimised away by the compiler
            if ωT <: AbstractMatrix
                ωi = pinv.ω
            else
                pnti = view(points, i, :)
                ωi = pinv.ω(pnti, t, p)
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
