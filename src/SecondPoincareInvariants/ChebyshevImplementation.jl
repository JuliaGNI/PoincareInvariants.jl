"""
    ChebyshevImplementation

implementation of computation of second Poincare invariant by approximating surface with
Chebyshev polynomials
"""
module ChebyshevImplementation

using ...PoincareInvariants: @argcheck

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
	phasecoeffs::AbstractArray, ∂xcoeffs::AbstractArray, ∂ycoeffs::AbstractArray
)
	∂xvals = invpaduatransform!(plan.∂xvals, plan.invpaduaplan, ∂xcoeffs)
	∂yvals = invpaduatransform!(plan.∂yvals, plan.invpaduaplan, ∂ycoeffs)

	@views for i in axes(plan.intvals, 1)
		plan.intvals[i] = dot(∂yvals[i, :], Ω(phasepoints[i, :], t, p), ∂xvals[i, :])
	end

	paduatransform!(intcoeffs, plan.paduaplan, plan.intvals)

	intcoeffs
end

# struct ChebyshevPlan{T, I, P<:PaduaTransformPlan, DO, IO}
#     degree::Int
#     paduaplan::PP
#     phasecoeffs::Array{T, 3}
#     diffplan::DP
#     ∂x::Array{T, 3}
#     ∂y::Array{T, 3}
#     intplan::IP  # getting coefficients of integrand to integrate
#     intcoeffs::Matrix{T}
#     integrator::Vector{T}
# end

# function _compute!(P::ChebyshevPlan, Ω, D, phasepoints)
#     paduatransform!(P.phasecoeffs, P.paduaplan, phasepoints)
#     differentiate!(P.∂x, P.∂y, P.diffplan, P.phasecoeffs)
#     getintegrand!(P.intcoeffs, P.intcoeffsplan, Ω, phasepoints, P.phasecoeffs, P.∂x, P.∂y)
#     integrate(P.intcoeffs, P.intplan)
# end

end  # module ChebyshevImplementation
