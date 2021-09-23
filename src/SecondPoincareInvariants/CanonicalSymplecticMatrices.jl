@reexport module CanonicalSymplecticMatrices

import Base

export CanonicalSymplecticMatrix

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

end  # module
