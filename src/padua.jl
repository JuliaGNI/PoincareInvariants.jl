# For copy of paduatransform!
using LinearAlgebra: rmul!
using FastTransforms: paduavalsmat, trianglecfsvec!

"""
    get_padua_points([T::Type{<:Real}, ]N::Integer)::Matrix{T}

returns `N` padua points on square `0..1 × 0..1` as matrix with elment type `T`. Each row
represens a point.
"""
function get_padua_points(::Type{T}, padua_num::Integer)::Matrix{T} where T <: Real
    check_padua_num(padua_num)
    out = Matrix{T}(undef, padua_num, 2)
    n = Int(get_n(padua_num))
    m = 0
    delta = 0
    NN = fld(n + 2, 2)
    @inbounds for k = n:-1:0
        if isodd(n)
            delta = mod(k, 2)
        end
        @inbounds for j = NN+delta:-1:1
            m += 1

            out[m, 1] = (sinpi(T(k) / T(n) - T(0.5)) + 1) / 2

            a = isodd(n - k) ? 1 : 2
            out[m, 2] = (sinpi(T(2j - a) / T(n + 1) - T(0.5)) + 1) / 2
        end
    end
    return out
end

get_padua_points(padua_num::Integer) = get_padua_points(Float64, padua_num)

function padua_transform!(
    out::AbstractVector{T},
    v::AbstractVector{T},
    P::PaduaTransformPlan
) where T
    axes(out, 1) == axes(v, 1) || error("axes of input and output vectors must match")
    N=length(v)
    n=Int(cld(-3+sqrt(1+8N),2))
    vals=paduavalsmat(P,v)
    tensorcfs = P.dctplan*vals
    m,l=size(tensorcfs)
    rmul!(tensorcfs,T(2)/(n*(n+1)))
    rmul!(view(tensorcfs,1,:),0.5)
    rmul!(view(tensorcfs,:,1),0.5)
    rmul!(view(tensorcfs,m,:),0.5)
    rmul!(view(tensorcfs,:,l),0.5)
    trianglecfsvec!(out,P,tensorcfs)
end
