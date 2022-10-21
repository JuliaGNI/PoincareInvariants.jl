module FirstFinDiffPlans

using LinearAlgebra: dot
using StaticArrays: SVector

using ..PoincareInvariants: FirstPoincareInvariant, getpointnum, getform

import ..PoincareInvariants: compute!, getpoints, getpointspec

struct FirstFinDiffPlan{T, D} end

FirstFinDiffPlan{T, D}(θ, N) where {T, D} = FirstFinDiffPlan{T, D}()

function compute!(
    pinv::FirstPoincareInvariant{T, D, <:Any, <:FirstFinDiffPlan}, t::Real, p
) where {T, D}
    zs = pinv.points
    N = getpointnum(pinv)
    θ = getform(pinv)
    I = zero(T)

    dzi = Vector{T}(undef, D)

    for i in 2:N-1
        zi = view(zs, i, :)
        dzi .= @views (zs[i+1, :] .- zs[i-1, :]) ./ 2
        I += dot(θ(zi, t, p), dzi)
    end

    # special cases i = 1
    z1 = view(zs, 1, :)
    dzi .= @views (zs[2, :] .- zs[N, :]) ./ 2
    I += dot(θ(z1, t, p), dzi)

    # special cases i = N
    zN = view(zs, N, :)
    dzi .= @views (zs[1, :] .- zs[N-1, :]) ./ 2
    I += dot(θ(zN, t, p), dzi)

    return I
end

## getpoints and getpointspec ##

# must be even to work with simpson weights
function getpointspec(N::Integer, ::Type{<:FirstFinDiffPlan})::Int
    if N < 2
        return 2
    elseif isodd(N)
        return N + 1
    else
        return N
    end
end

function getpoints(f, ::Type{T}, N::Integer, ::Type{<:FirstFinDiffPlan}) where T
    D = length(f(zero(T)))
    out = Matrix{T}(undef, N, D)

    for (i, x) in enumerate(range(0, 1, length=N+1)[1:end-1])
        out[i, :] .= f(x)
    end

    # Should return vector instead of matrix in 1D case
    return D == 1 ? vec(out) : out
end

end  # module FirstFinDiffPlans
