"""
    SecondChebyshevPlans

implementation of computation of second Poincare invariant by approximating surface with
Chebyshev polynomials
"""
module SecondChebyshevPlans

import ..PoincareInvariants: compute!, getpoints, getpointspec

using ..PoincareInvariants: SecondPoincareInvariant, @argcheck

using LinearAlgebra

using ChebyshevTransforms

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

## getintegrand for ω Callable and out-of-place ##

struct CallIntPlan{T, D, IP, P}
    invpaduaplan::IP
    ∂xvals::Matrix{T}
    ∂yvals::Matrix{T}
    intvals::Vector{T}
    paduaplan::P
end

function CallIntPlan{T, D}(degree) where {T, D}
    invpaduaplan = InvPaduaTransformPlan{Float64}(degree)
    ∂xvals = Matrix{T}(undef, getpaduanum(degree), D)
    ∂yvals = Matrix{T}(undef, getpaduanum(degree), D)
    intvals = Vector{T}(undef, getpaduanum(degree))
    paduaplan = PaduaTransformPlan{Float64}(degree)

    CallIntPlan{T, D, typeof(invpaduaplan), typeof(paduaplan)}(
        invpaduaplan, ∂xvals, ∂yvals, intvals, paduaplan
    )
end

function getintegrand!(
    intcoeffs::AbstractMatrix, plan::CallIntPlan{T, D}, ω::ωT,
    points, t, p, ∂xcoeffs, ∂ycoeffs
) where {T, D, ωT}
    invpaduatransform!(plan.∂xvals, plan.invpaduaplan, ∂xcoeffs)
    invpaduatransform!(plan.∂yvals, plan.invpaduaplan, ∂ycoeffs)

    for i in axes(plan.intvals, 1)
        # This if statement should hopefully get optimised away by the compiler
        if ωT <: AbstractMatrix
            ωi = ω
        else
            pnti = view(points, i, :)
            ωi = ω(pnti, t, p)
        end

        ∂xi = view(plan.∂xvals, i, :)
        ∂yi = view(plan.∂yvals, i, :)
        plan.intvals[i] = dot(∂yi, ωi, ∂xi)
    end

    paduatransform!(intcoeffs, plan.paduaplan, plan.intvals)

    intcoeffs
end

## SecondChebyshevPlan and compute! ##

getintplan(::Type{T}, ::Any, ::Val{D}, degree) where {T, D} = CallIntPlan{T, D}(degree)
# getintplan(::Type{T}, ω::AbstractMatrix, D, degree) where T = ConstIntPlan{T}(ω, D, degree)

struct SecondChebyshevPlan{T, D, IP, PP<:PaduaTransformPlan}
    degree::Int
    paduaplan::PP
    phasecoeffs::NTuple{D, Matrix{T}}
    diffplan::DiffPlan{T}
    ∂x::NTuple{D, Matrix{T}}
    ∂y::NTuple{D, Matrix{T}}
    intplan::IP  # getting coefficients of integrand to integrate
    intcoeffs::Matrix{T}
    intweights::Vector{T}
end

function SecondChebyshevPlan{T, D}(ω, N::Integer) where {T, D}
    degree = nextdegree(N)

    paduaplan = PaduaTransformPlan{Float64}(degree)
    phasecoeffs = ntuple(_ -> zeros(T, degree+1, degree+1), D)
    diffplan = DiffPlan{T}(degree)
    ∂x = ntuple(_ -> Matrix{T}(undef, degree+1, degree+1), D)
    ∂y = ntuple(_ -> Matrix{T}(undef, degree+1, degree+1), D)
    intplan = getintplan(T, ω, Val(D), degree)
    intcoeffs = zeros(T, degree+1, degree+1)
    intweights = getintweights(T, degree)

    SecondChebyshevPlan{T, D, typeof(intplan), typeof(paduaplan)}(
        degree, paduaplan, phasecoeffs, diffplan, ∂x, ∂y, intplan, intcoeffs, intweights
    )
end

function compute!(
    pinv::SecondPoincareInvariant{<:Any, <:Any, <:Any, <:Any, P}, t::Real, p
) where P <: SecondChebyshevPlan
    plan = pinv.plan
    points = pinv.points
    paduatransform!(plan.phasecoeffs, plan.paduaplan, points)
    differentiate!(plan.∂x, plan.∂y, plan.diffplan, plan.phasecoeffs)
    getintegrand!(plan.intcoeffs, plan.intplan, pinv.ω, points, t, p, plan.∂x, plan.∂y)
    integrate(plan.intcoeffs, plan.intweights)
end

# function _compute!(plan::SecondChebyshevPlan, ω::AbstractMatrix, D, points)
#     paduatransform!(plan.phasecoeffs, plan.paduaplan, points)
#     differentiate!(plan.∂x, plan.∂y, plan.diffplan, plan.phasecoeffs)
#     getintegrand!(plan.intcoeffs, plan.intplan, ω, plan.∂x, plan.∂y)
#     integrate(plan.intcoeffs, plan.intplan)
# end

## getpoints, getpointspec and getpointnum ##

getpointspec(N::Integer, ::Type{<:SecondChebyshevPlan}) = nextpaduanum(N)

function getpoints(f, ::Type{T}, N, ::Type{<:SecondChebyshevPlan}) where T
    degree = getdegree(getpointspec(N, SecondChebyshevPlan))
    return getpaduapoints(T, degree) do x, y
        f((x + T(1)) / T(2), (y + T(1)) / T(2))
    end
end

end  # module Chebyshev
