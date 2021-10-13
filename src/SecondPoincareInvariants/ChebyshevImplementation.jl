"""
    ChebyshevImplementation

implementation of computation of second Poincare invariant by approximating surface with
Chebyshev polynomials
"""
module ChebyshevImplementation

using ...PoincareInvariants: @argcheck
import ...PoincareInvariants: compute!, getpoints, getpointnum

using Base: Callable
using LinearAlgebra

include("PaduaTransforms.jl")
using .PaduaTransforms

const AbstractArray3 = AbstractArray{<:Any, 3}

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

function differentiate!(∂x::AbstractArray3, ∂y::AbstractArray3, P::DiffPlan, coeffs::AbstractArray3)
    axes(∂x, 3) == axes(∂y, 3) == axes(coeffs, 3) || error()

    @views for i in axes(coeffs, 3)
        differentiate!(∂x[:, :, i], ∂y[:, :, i], P, coeffs[:, :, i])
    end

    ∂x, ∂y
end

## Integration ##

getintegrator(::Type{T}, n) where T = T[isodd(i) ? 0 : T(2) / T(1 - i^2) for i in 0:n]
getintegrator(n) = getintegrator(Float64, n)

integrate(coeffs, integrator) = dot(integrator, coeffs, integrator)

## getintegrand for Ω Callable and out-of-place ##

struct OOPIntPlan{T, IP, P}
    invpaduaplan::IP
    ∂xvals::Matrix{T}
    ∂yvals::Matrix{T}
    intvals::Vector{T}
    paduaplan::P
end

function OOPIntPlan{T}(D, degree) where T
    invpaduaplan = InvPaduaTransformPlan{T}(degree)
    ∂xvals = Matrix{T}(undef, getpaduanum(degree), D)
    ∂yvals = Matrix{T}(undef, getpaduanum(degree), D)
    intvals = Vector{T}(undef, getpaduanum(degree))
    paduaplan = PaduaTransformPlan{T}(degree)

    OOPIntPlan{T, typeof(invpaduaplan), typeof(paduaplan)}(
        invpaduaplan, ∂xvals, ∂yvals, intvals, paduaplan
    )
end

function getintegrand!(
    intcoeffs::AbstractMatrix, plan::OOPIntPlan, Ω::Callable, phasepoints::AbstractMatrix, t, p,
    ∂xcoeffs::AbstractArray3, ∂ycoeffs::AbstractArray3
)
    ∂xvals = invpaduatransform!(plan.∂xvals, plan.invpaduaplan, ∂xcoeffs)
    ∂yvals = invpaduatransform!(plan.∂yvals, plan.invpaduaplan, ∂ycoeffs)

    @views for i in axes(plan.intvals, 1)
        plan.intvals[i] = dot(∂yvals[i, :], Ω(phasepoints[i, :], t, p), ∂xvals[i, :])
    end

    paduatransform!(intcoeffs, plan.paduaplan, plan.intvals)

    intcoeffs
end

## ChebyshevPlan and _compute! ##

getintplan(::Type{T}, ::Callable, D, degree, ::Val{false}) where T = OOPIntPlan{T}(D, degree)
# getintplan(::Type{T}, ::Callable, D, degree, ::Val{true}) where T = IPIntPlan{T}(D, degree)
# getintplan(::Type{T}, Ω::AbstractMatrix, D, degree, ::Val{nothing}) where T = IPIntPlan{T}(Ω, D, degree)

struct ChebyshevPlan{T, IP, PP<:PaduaTransformPlan}
    degree::Int
    paduaplan::PP
    phasecoeffs::Array{T, 3}
    diffplan::DiffPlan{T}
    ∂x::Array{T, 3}
    ∂y::Array{T, 3}
    intplan::IP  # getting coefficients of integrand to integrate
    intcoeffs::Matrix{T}
    integrator::Vector{T}
end

function ChebyshevPlan{T}(Ω::Callable, D::Integer, N::Integer, ::Val{inplace}) where {T, inplace}
    degree = getdegree(nextpaduanum(N))

    paduaplan = PaduaTransformPlan{T}(degree)
    phasecoeffs = zeros(T, degree+1, degree+1, D)
    diffplan = DiffPlan{T}(degree)
    ∂x = Array{T, 3}(undef, degree+1, degree+1, D)
    ∂y = Array{T, 3}(undef, degree+1, degree+1, D)
    intplan = getintplan(T, Ω, D, degree, Val(inplace))
    intcoeffs = zeros(T, degree+1, degree+1)
    integrator = getintegrator(T, degree)

    ChebyshevPlan{T, typeof(intplan), typeof(paduaplan)}(
        degree, paduaplan, phasecoeffs, diffplan, ∂x, ∂y, intplan, intcoeffs, integrator
    )
end

function compute!(plan::ChebyshevPlan, Ω::Callable, phasepoints::AbstractMatrix, t, p)
    paduatransform!(plan.phasecoeffs, plan.paduaplan, phasepoints)
    differentiate!(plan.∂x, plan.∂y, plan.diffplan, plan.phasecoeffs)
    getintegrand!(plan.intcoeffs, plan.intplan, Ω, phasepoints, t, p, plan.∂x, plan.∂y)
    integrate(plan.intcoeffs, plan.integrator)
end

# function _compute!(plan::ChebyshevPlan, Ω::AbstractMatrix, D, phasepoints)
#     paduatransform!(plan.phasecoeffs, plan.paduaplan, phasepoints)
#     differentiate!(plan.∂x, plan.∂y, plan.diffplan, plan.phasecoeffs)
#     getintegrand!(plan.intcoeffs, plan.intplan, Ω, plan.∂x, plan.∂y)
#     integrate(plan.intcoeffs, plan.intplan)
# end

## getpoints and getpointnum ##

getpointnum(plan::ChebyshevPlan) = getpaduanum(plan.degree)
getpoints(plan::ChebyshevPlan) = getpaduapoints(plan.degree) do x, y
    (x + 1) / 2, (y + 1) / 2
end

getpoints(f::Function, plan::ChebyshevPlan) = getpaduapoints(plan.degree) do x, y
    f((x + 1) / 2, (y + 1) / 2)
end

end  # module ChebyshevImplementation
