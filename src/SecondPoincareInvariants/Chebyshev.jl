"""
    Chebyshev

implementation of computation of second Poincare invariant by approximating surface with
Chebyshev polynomials
"""
module Chebyshev

import ...PoincareInvariants: compute!, getpoints, getpointnum
import ..SecondPoincareInvariants: getpointspec

using ...PoincareInvariants: @argcheck
using ..SecondPoincareInvariants: SecondPoincareInvariant

using Base: Callable
using LinearAlgebra

include("PaduaTransforms.jl")
using .PaduaTransforms

## Differentiation ##

function getdiffmat(::Type{T}, degree::Integer) where T
    D = zeros(T, degree+1, degree+1)

    D[1, 2:2:end] .= 1:2:degree

    for i in 3:2:degree+1
        D[2:2:i-1, i] .= 2 * (i-1)
    end

    for i in 4:2:degree+1
        D[3:2:i-1, i] .= 2 * (i-1)
    end

    D
end

struct DiffPlan{T}
    D::Matrix{T}
end

DiffPlan{T}(degree::Integer) where T = DiffPlan{T}(getdiffmat(T, degree))

# differentiate in the Chebyshev basis
function differentiate!(∂x::AbstractMatrix, ∂y::AbstractMatrix, P::DiffPlan, coeffs::AbstractMatrix)
    ∂x[:, :] = coeffs  # differentiate each row
    rmul!(∂x, LowerTriangular(P.D'))

    ∂y[:, :] = coeffs  # differentiate each column
    lmul!(UpperTriangular(P.D), ∂y)

    ∂x, ∂y
end

function differentiate!(∂x, ∂y, P::DiffPlan, coeffs)
    eltype(∂x) <: AbstractMatrix && eltype(∂y) <: AbstractMatrix &&
        eltype(coeffs) <: AbstractMatrix || throw(ArgumentError(
            "coefficients must be AbstractMatrix or iterable thereof"))
    length(∂x) == length(∂y) == length(coeffs) || throw(ArgumentError(
        "number of coefficient matrices must match number of derivative matrices to write to"))

    for (∂xi, ∂yi, coeffsi) in zip(∂x, ∂y, coeffs)
        differentiate!(∂xi, ∂yi, P, coeffsi)
    end

    ∂x, ∂y
end

## Integration ##

getintweights(::Type{T}, n) where T = T[isodd(i) ? 0 : T(2) / T(1 - i^2) for i in 0:n]
getintweights(n) = getintweights(Float64, n)

integrate(coeffs, intweights) = dot(intweights, coeffs, intweights)

## getintegrand for Ω Callable and out-of-place ##

struct CallIntPlan{T, IP, P}
    invpaduaplan::IP
    phasevals::Matrix{T}
    ∂xvals::Matrix{T}
    ∂yvals::Matrix{T}
    intvals::Vector{T}
    paduaplan::P
end

function CallIntPlan{T}(D, degree) where T
    invpaduaplan = InvPaduaTransformPlan{T}(degree)
    phasevals = Matrix{T}(undef, D, getpaduanum(degree))
    ∂xvals = Matrix{T}(undef, D, getpaduanum(degree))
    ∂yvals = Matrix{T}(undef, D, getpaduanum(degree))
    intvals = Vector{T}(undef, getpaduanum(degree))
    paduaplan = PaduaTransformPlan{T}(degree)

    CallIntPlan{T, typeof(invpaduaplan), typeof(paduaplan)}(
        invpaduaplan, phasevals, ∂xvals, ∂yvals, intvals, paduaplan
    )
end

function getintegrand!(
    intcoeffs::AbstractMatrix, plan::CallIntPlan{T}, Ω::Callable,
    phasepoints, t, p, ∂xcoeffs, ∂ycoeffs
) where T
    invpaduatransform!(eachrow(plan.∂xvals), plan.invpaduaplan, ∂xcoeffs)
    invpaduatransform!(eachrow(plan.∂yvals), plan.invpaduaplan, ∂ycoeffs)

    D = length(phasepoints)
    for d in 1:D
        plan.phasevals[d, :] .= phasepoints[d]
    end

    for i in axes(plan.intvals, 1)
        pnti = view(plan.phasevals, :, i)
        ∂xi = view(plan.∂xvals, :, i)
        ∂yi = view(plan.∂yvals, :, i)
        plan.intvals[i] = dot(∂yi, Ω(pnti, t, p), ∂xi)
    end

    paduatransform!(intcoeffs, plan.paduaplan, plan.intvals)

    intcoeffs
end

## ChebyshevPlan and compute! ##

getintplan(::Type{T}, ::Callable, D, degree) where T = CallIntPlan{T}(D, degree)
# getintplan(::Type{T}, Ω::AbstractMatrix, D, degree) where T = ConstIntPlan{T}(Ω, D, degree)

struct ChebyshevPlan{T, IP, PP<:PaduaTransformPlan}
    degree::Int
    paduaplan::PP
    phasecoeffs::Vector{Matrix{T}}
    diffplan::DiffPlan{T}
    ∂x::Vector{Matrix{T}}
    ∂y::Vector{Matrix{T}}
    intplan::IP  # getting coefficients of integrand to integrate
    intcoeffs::Matrix{T}
    intweights::Vector{T}
end

function ChebyshevPlan{T}(Ω::Callable, D::Integer, N::Integer) where T
    degree = getdegree(nextpaduanum(N))

    paduaplan = PaduaTransformPlan{T}(degree)
    phasecoeffs = [zeros(T, degree+1, degree+1) for _ in 1:D]
    diffplan = DiffPlan{T}(degree)
    ∂x = [Matrix{T}(undef, degree+1, degree+1) for _ in 1:D]
    ∂y = [Matrix{T}(undef, degree+1, degree+1) for _ in 1:D]
    intplan = getintplan(T, Ω, D, degree)
    intcoeffs = zeros(T, degree+1, degree+1)
    intweights = getintweights(T, degree)

    ChebyshevPlan{T, typeof(intplan), typeof(paduaplan)}(
        degree, paduaplan, phasecoeffs, diffplan, ∂x, ∂y, intplan, intcoeffs, intweights
    )
end

function compute!(
    pinv::SecondPoincareInvariant{T, ΩT, NT, P}, phasepoints, t, p
) where {T, ΩT <: Callable, NT, P <: ChebyshevPlan}
    plan = pinv.plan
    paduatransform!(plan.phasecoeffs, plan.paduaplan, phasepoints)
    differentiate!(plan.∂x, plan.∂y, plan.diffplan, plan.phasecoeffs)
    getintegrand!(plan.intcoeffs, plan.intplan, pinv.Ω, phasepoints, t, p, plan.∂x, plan.∂y)
    integrate(plan.intcoeffs, plan.intweights)
end

# function _compute!(plan::ChebyshevPlan, Ω::AbstractMatrix, D, phasepoints)
#     paduatransform!(plan.phasecoeffs, plan.paduaplan, phasepoints)
#     differentiate!(plan.∂x, plan.∂y, plan.diffplan, plan.phasecoeffs)
#     getintegrand!(plan.intcoeffs, plan.intplan, Ω, plan.∂x, plan.∂y)
#     integrate(plan.intcoeffs, plan.intplan)
# end

## getpoints, getpointspec and getpointnum ##

getpointnum(N::Integer, ::Type{<:ChebyshevPlan}) = nextpaduanum(N)
getpointnum(dims::NTuple{2, Integer}, ::Type{<:ChebyshevPlan}) = nextpaduanum(dims[1] * dims[2])

function getpoints(f, ::Type{T}, N, ::Type{<:ChebyshevPlan}) where T
    degree = nextdegree(N)
    return getpaduapoints(T, degree) do x, y
        f((x + T(1)) / T(2), (y + T(1)) / T(2))
    end
end

end  # module Chebyshev
