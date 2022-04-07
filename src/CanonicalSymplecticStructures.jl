module CanonicalSymplecticStructures

import Base
import LinearAlgebra

export CanonicalSymplecticMatrix

"""
    CanonicalSymplecticMatrix{T}(n::Integer)

constructs a canonical symplectic matrix of size `(n, n)` with eltype `T`.
`n` must be even. See the examples to see the form of the canonical symplectic matrix as
defined here.

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
        iseven(n) || throw(ArgumentError("CanonicalSymplecticMatrix must have even size"))
        new{T}(n ÷ 2)
    end
end

CanonicalSymplecticMatrix(n::Integer) = CanonicalSymplecticMatrix{Int}(n)

Base.size(C::CanonicalSymplecticMatrix) = (sz = C.mid * 2; (sz, sz))

function Base.getindex(C::CanonicalSymplecticMatrix{T}, i1::Int, i2::Int) where T
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
    axes(x, 1) == axes(y, 1) == 1:2m || error()
    s = zero(promote_type(eltype(x), eltype(y)))
    @inbounds for i in 1:m
        s += y[i] * x[m + i]
        s -= y[m + i] * x[i]
    end

    s
end

end  # module
