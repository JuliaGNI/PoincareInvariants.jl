get_n(padua_num) = (sqrt(1 + 8padua_num) - 3) / 2

get_padua_num(n) = (n + 1) * (n + 2) ÷ 2

function check_padua_num(padua_num)
    if !isinteger(get_n(padua_num))
        throw(ArgumentError("number of padua points must equal (n + 1) * (n + 2) ÷ 2"))
    end
end

"""
    get_padua_points(N::Integer)::AbstractMatrix{Float64}

returns `N` padua points on square `0..1 × 0..1`.
"""
function get_padua_points(padua_num::Integer)::Vector{SVector{2, Float64}}
    check_padua_num(padua_num)
    out = Vector{SVector{2, Float64}}(undef, padua_num)
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

            v1 = (sinpi(k / n - 0.5) + 1) / 2

            a = isodd(n - k) ? 1 : 2
            v2 = (sinpi((2j - a) / (n + 1) - 0.5) + 1) / 2

            out[m] = SVector{2, Float64}(v1, v2)
        end
    end
    return out
end

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
