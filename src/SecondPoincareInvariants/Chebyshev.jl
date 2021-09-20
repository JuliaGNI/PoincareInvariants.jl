# For copy of paduatransform!
using LinearAlgebra: rmul!
using FastTransforms: paduavalsmat, trianglecfsvec!

getdegree(coeffnum) = (sqrt(1 + 8coeffnum) - 3) / 2
getcoeffnum(degree) = (degree + 1) * (degree + 2) ÷ 2
getpointnum(degree) = (degree + 1)^2

getpaduanum(degree) = getcoeffnum(degree)
nextpaduanum(N) = getpaduanum(ceil(Int, getdegree(N)))

function checkpaduanum(paduanum)
    if !isinteger(getdegree(paduanum))
        throw(ArgumentError("number of Padua points / coeffs must be a triangle number 1, 3, 6, 10, 15..."))
    end
end

"""
    getpaduapoints([T::Type{<:Real}, ]n::Integer)::Matrix{T}

returns padua points corresponding to degree `n` chebyshev polynomial on square `0..1 × 0..1`
as matrix with element type `T`. Each row represens a point.
"""
function getpaduapoints(::Type{T}, n::Integer)::Matrix{T} where T <: Real
    paduanum = getpaduanum(n)
    out = Matrix{T}(undef, paduanum, 2)
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

getpaduapoints(n::Integer) = getpaduapoints(Float64, n)

function paduatransform!(
    out::AbstractVector{T},
    v::AbstractVector{T},
    P::PaduaTransformPlan
) where T
    axes(out, 1) == axes(v, 1) || error("axes of input and output vectors must match")
    N = length(v)
    checkpaduanum(N)
    n = Int(getdegree(N))
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
