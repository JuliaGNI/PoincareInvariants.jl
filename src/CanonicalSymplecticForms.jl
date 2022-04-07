module CanonicalSymplecticForms

import Base
import LinearAlgebra

export CanonicalSymplecticTwoForm

"""
    CanonicalSymplecticTwoForm{T}(n::Integer)

constructs a canonical symplectic matrix of size `(n, n)` with eltype `T`.
`n` must be even. See the examples to see the form of the canonical symplectic matrix as
defined here.

# Examples

```jldoctest
julia> CanonicalSymplecticTwoForm(4)
4×4 CanonicalSymplecticTwoForm{Int64}:
 0  0  -1   0
 0  0   0  -1
 1  0   0   0
 0  1   0   0

julia> CanonicalSymplecticTwoForm{Int32}(6)
6×6 CanonicalSymplecticTwoForm{Int32}:
 0  0  0  -1   0   0
 0  0  0   0  -1   0
 0  0  0   0   0  -1
 1  0  0   0   0   0
 0  1  0   0   0   0
 0  0  1   0   0   0
```
"""
struct CanonicalSymplecticTwoForm{T} <: AbstractMatrix{T}
    mid::Int
    function CanonicalSymplecticTwoForm{T}(n::Integer) where T
        iseven(n) || throw(ArgumentError("CanonicalSymplecticTwoForm must have even size"))
        new{T}(n ÷ 2)
    end
end

CanonicalSymplecticTwoForm(n::Integer) = CanonicalSymplecticTwoForm{Int}(n)

Base.size(C::CanonicalSymplecticTwoForm) = (sz = C.mid * 2; (sz, sz))

function Base.getindex(C::CanonicalSymplecticTwoForm{T}, i1::Int, i2::Int) where T
    if i2 - C.mid == i1
        return T(-1)
    elseif i1 - C.mid == i2
        return T(1)
    else
        return T(0)
    end
end

function LinearAlgebra.dot(x::AbstractVector, C::CanonicalSymplecticTwoForm, y::AbstractVector)
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
