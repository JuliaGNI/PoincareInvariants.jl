module FirstFourierPlans

using FFTW
using LinearAlgebra: mul!

using ..PoincareInvariants: FirstPoincareInvariant, getpointnum, getform, getplan

import ..PoincareInvariants: compute!, getpoints, getpointspec

export FirstFourierPlan

struct FirstFourierPlan{T, D, FTP}
    θs::Matrix{T}
    input::Vector{Float64}
    θks::Vector{Complex{Float64}}
    zks::Vector{Complex{Float64}}
    fftplan::FTP
end

function FirstFourierPlan{T, D}(θ, N) where {T, D}
    θs = Matrix{T}(undef, N, D)
    input = Vector{Float64}(undef, N)
    θks = Vector{Complex{Float64}}(undef, N ÷ 2 + 1)
    zks = Vector{Complex{Float64}}(undef, N ÷ 2 + 1)

    # this makes use of FFTW internals, so it may easily break
    plan = FFTW.rFFTWPlan{Float64, FFTW.FORWARD, false, 1}(
        input, θks, 1:1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT
    )
    FirstFourierPlan{T, D, typeof(plan)}(θs, input, θks, zks, plan)
end

function compute!(
    pinv::FirstPoincareInvariant{T, D, <:Any, <:FirstFourierPlan}, t::Real, p
) where {T, D}

    zs = pinv.points
    N = getpointnum(pinv)
    plan = getplan(pinv)

    θ = getform(pinv)
    for (i, z) in enumerate(eachrow(zs))
        plan.θs[i, :] .= θ(z, t, p)
    end

    I = zero(T)

    for d in 1:D
        copyto!(plan.input, view(plan.θs, :, d))
        mul!(plan.θks, plan.fftplan, plan.input)

        copyto!(plan.input, view(zs, :, d))
        mul!(plan.zks, plan.fftplan, plan.input)

        for k in 1:(N ÷ 2 + 1)
            I += T(real(plan.zks[k] * conj(plan.θks[k]) * im * 4π * (k-1) / N^2))
        end
    end

    return I
end

## getpoints and getpointspec ##

function getpointspec(N::Integer, ::Type{<:FirstFourierPlan})::Int
    N < 0 && error()
    return N
end

function getpoints(f, ::Type{T}, N::Integer, ::Type{<:FirstFourierPlan}) where T
    D = length(f(zero(T)))
    out = Matrix{T}(undef, N, D)

    for (i, x) in enumerate(range(0, 1, length=N+1)[1:end-1])
        out[i, :] .= f(x)
    end

    # Should return vector instead of matrix in 1D case
    return D == 1 ? vec(out) : out
end

end  # module FirstFourierPlans
