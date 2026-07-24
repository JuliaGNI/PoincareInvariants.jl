module CanonicalSymplecticForms

import Base
import LinearAlgebra

export canonical_one_form!, CanonicalSymplecticMatrix

function checkn(T, n::Integer)
    n > 0 || throw(ArgumentError("$T must have positive size"))
    iseven(n) || throw(ArgumentError("$T must have even size"))
end

"""
    canonical_one_form!(out, t, z, p)

writes the canonical one form ``\\vartheta`` at the phase space point `z` into `out`, following
the in-place `form(out, t, z, p)` convention. For a phase space point `z = (q, p)` of even
length `n = 2m`, the canonical one form is ``\\vartheta = (p, 0)``, i.e. `out[1:m] .= z[m+1:n]`
and `out[m+1:n] .= 0`. The time `t` and parameters `p` are ignored.
"""
function canonical_one_form!(out, t, z, p)
    n = length(z)
    iseven(n) || throw(ArgumentError("z must have even length"))
    mid = n ÷ 2
    @inbounds @views begin
        out[1:mid] .= z[mid+1:n]
        out[mid+1:n] .= zero(eltype(out))
    end
    return out
end

"""
    CanonicalSymplecticMatrix{T}(n::Integer)

constructs a canonical symplectic matrix of size `(n, n)` with eltype `T`.
`n` must be even and positive. See the examples to see the form of the canonical symplectic
matrix as defined here.

# Examples

```jldoctest
julia> CanonicalSymplecticMatrix(4)
4×4 CanonicalSymplecticMatrix{Int64}:
 0  0  -1   0
 0  0   0  -1
 1  0   0   0
 0  1   0   0

julia> CanonicalSymplecticMatrix{Int32}(6)
6×6 CanonicalSymplecticMatrix{Int32}:
 0  0  0  -1   0   0
 0  0  0   0  -1   0
 0  0  0   0   0  -1
 1  0  0   0   0   0
 0  1  0   0   0   0
 0  0  1   0   0   0
```
"""
struct CanonicalSymplecticMatrix{T} <: AbstractMatrix{T}
    mid::Int
    function CanonicalSymplecticMatrix{T}(n::Integer) where T
        checkn(CanonicalSymplecticMatrix, n)
        new{T}(n ÷ 2)
    end
end

CanonicalSymplecticMatrix(n::Integer) = CanonicalSymplecticMatrix{Int}(n)

Base.size(C::CanonicalSymplecticMatrix) = (sz = C.mid * 2; (sz, sz))

function Base.getindex(C::CanonicalSymplecticMatrix{T}, i1::Int, i2::Int) where T
    @boundscheck let n = C.mid * 2
        if !(1 ≤ i1 ≤ n && 1 ≤ i2 ≤ n)
            msg = "attempt to access $n-element CanonicalSymplecticMatrix{$T} at index [$i1, $i2]"
            throw(BoundsError(msg))
        end
    end

    if i2 - C.mid == i1
        return T(-1)
    elseif i1 - C.mid == i2
        return T(1)
    else
        return T(0)
    end
end

function LinearAlgebra.dot(x::AbstractVector, C::CanonicalSymplecticMatrix, y::AbstractVector)
    m = C.mid
    axes(x, 1) == axes(y, 1) == 1:2m || throw(DimensionMismatch())
    s = zero(promote_type(eltype(x), eltype(y)))
    @inbounds for i in 1:m
        s += y[i] * x[m + i]
        s -= y[m + i] * x[i]
    end

    s
end

end  # module
