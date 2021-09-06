_to_n(N) = (-3 + sqrt(1 + 8N)) / 2
_to_N(n) = (n + 1) * (n + 2) ÷ 2

# get minimum allowed padua point number ≥ N
_get_min_padua_num(N) = _to_N(ceil(Int, _to_n(N)))

function _padua_points(N::Integer)::Vector{SVector{2, Float64}}
    out = Vector{SVector{2, Float64}}(undef, N)
    n = Int(_to_n(N))
    m = 0
    delta = 0
    NN = fld(n + 2, 2)
    @inbounds for k = n:-1:0
        if isodd(n)
            delta = mod(k, 2)
        end
        @inbounds for j = NN+delta:-1:1
            m += 1

            v1 = sinpi(k / n - 0.5)

            a = isodd(n - k) ? 1 : 2
            v2 = sinpi((2j - a) / (n + 1) - 0.5)

            out[m] = SVector{2, Float64}(v1, v2)
        end
    end
    return out
end

"""
    get_padua_points(N::Integer)::Vector{SVector{2, Float64}}

returns minimum number of padua points ≥ `N` on square `-1..1 × -1..1`.
"""
get_padua_points(N::Integer) where T = _padua_points(_get_min_padua_num(N))

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
