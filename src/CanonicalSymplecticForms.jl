module CanonicalSymplecticForms

import Base
import LinearAlgebra

export CanonicalSymplecticOneForm, CanonicalSymplecticMatrix, CanonicalSymplecticTwoForm

function checkn(T, n::Integer)
    n > 0 || throw(ArgumentError("$T must have positive size"))
    iseven(n) || throw(ArgumentError("$T must have even size"))
end

struct CanonicalSymplecticVector{T, VT <: AbstractVector{T}} <: AbstractVector{T}
    p::VT

    function CanonicalSymplecticVector{T, VT}(p) where {T, VT}
        Base.require_one_based_indexing(p)
        new{T, VT}(p)
    end
end

CanonicalSymplecticVector{T}(p::VT) where {T, VT} = CanonicalSymplecticVector{T, VT}(p)
CanonicalSymplecticVector(p::VT) where VT = CanonicalSymplecticVector{eltype(p), VT}(p)

Base.size(C::CanonicalSymplecticVector) = (length(C.p) * 2,)

function Base.getindex(C::CanonicalSymplecticVector{T}, i::Int) where T
    mid = length(C.p)
    @boundscheck let n = mid * 2
        if !(1 ≤ i ≤ n)
            msg = "attempt to access $n-element CanonicalSymplecticVector at index [$i]"
            throw(BoundsError(msg))
        end
    end

    if i ≤ mid
        return @inbounds C.p[i]
    else
        return zero(T)
    end
end

function LinearAlgebra.dot(C::CanonicalSymplecticVector, x::AbstractVector)
    m = length(C.p)
    axes(x, 1) == 1:2m || throw(DimensionMismatch())

    s = zero(promote_type(eltype(x), eltype(C)))
    @inbounds for i in 1:m
        s += C.p[i] * x[i]
    end

    s
end

# LinearAlgebra.dot(x::AbstractVector, C::CanonicalSymplecticVector) = LinearAlgebra.dot(C, x)

function CanonicalSymplecticOneForm(z, t, p)
    n = length(z)
    iseven(n) || throw(ArgumentError("z must have even length"))
    mid = n ÷ 2
    CanonicalSymplecticVector(view(z, mid+1:n))
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

CanonicalSymplecticTwoForm(z, t, p) = CanonicalSymplecticMatrix{eltype(z)}(length(z))

end  # module
