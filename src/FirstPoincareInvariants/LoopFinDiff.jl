module LoopFinDiff

using LinearAlgebra: dot
using Base: Callable
using StaticArrays: SVector

using ...PoincareInvariants: getpointnum, getform
using ..FirstPoincareInvariants: FirstPoincareInvariant

import ...PoincareInvariants: compute!, getpoints, getpointspec

struct LoopFinDiffPlan end

function compute!(
    pinv::FirstPoincareInvariant{T, D, <:Callable, <:LoopFinDiffPlan}, zs, t, p
) where {T, D}

    N = getpointnum(pinv)
    θ = getform(pinv)
    I = zero(T)

    for i in 2:N-1
        zi = SVector{D}(ntuple(d -> zs[d][i], D))
        dzi = ntuple(D) do d
            (zs[d][i+1] - zs[d][i-1])
        end |> SVector{D}

        # simpson weights with factor 2 absorbed in derivative
        w = isodd(i) ? T(1) / 3 : T(2) / 3
        I += w * dot(θ(zi, t, p), dzi)
    end

    # special cases i = 1
    z1 = SVector{D}(ntuple(d -> zs[d][1], D))
    dz1 = ntuple(D) do d
        (zs[d][2] - zs[d][N])
    end |> SVector{D}

    I += dot(θ(z1, t, p), dz1) * T(1) / 3

    # special cases i = N
    zN = SVector{D}(ntuple(d -> zs[d][N], D))
    dzN = ntuple(D) do d
        (zs[d][1] - zs[d][N-1])
    end |> SVector{D}

    I += dot(θ(zN, t, p), dzN) * T(2) / 3

    return I
end

## getpoints and getpointspec ##

function getpointspec(N::Integer, ::Type{<:LoopFinDiffPlan})::Int
    if N < 2
        return 2
    elseif isodd(N)
        return N + 1
    else
        return N
    end
end

function getpoints(f, ::Type{T}, N::Integer, ::Type{<:LoopFinDiffPlan}) where T
    D = length(f(zero(T)))
    out = ntuple(_ -> Vector{T}(undef, N), D)

    for (i, x) in enumerate(range(0, 1, length=N+1)[1:end-1])
        fpnt = f(x)
        for d in 1:D
            out[d][i] = fpnt[d]
        end
    end

    # Should return vector [...] instead of ([...],) in 1D case
    return D == 1 ? out[1] : out
end

end  # module LoopFinDiff
