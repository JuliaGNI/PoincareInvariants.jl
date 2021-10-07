"""
    ChebyshevImplementation

implementation of computation of second Poincare invariant by approximating surface with
Chebyshev polynomials
"""
module ChebyshevImplementation

using ...PoincareInvariants: @argcheck

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

# Differentiation in the Chebyshev basis
function differentiate!(∂x, ∂y, P, coeffs)
	# mul!(out, ::Matrix, ::LowerTriangular) is really slow,
	# so we use the standard method here
	mul!(∂x, coeffs, P.D')  # apply row by row
	mul!(∂y, UpperTriangular(P.D), coeffs)  # apply column by column

	∂x, ∂y
end

## Integration ##

getintegrator(::Type{T}, n) where T = T[isodd(i) ? 0 : T(2) / T(1 - i^2) for i in 0:n]
getintegrator(n) = getintegrator(Float64, n)

integrate(coeffs, integrator) = dot(integrator, coeffs, integrator)

# struct ChebyshevPlan{T, I, P<:PaduaTransformPlan, DO, IO}
#     degree::Int
#     paduaplan::PP
#     phasecoeffs::Array{T, 3}
#     diffplan::DP
#     ∂x::Array{T, 3}
#     ∂y::Array{T, 3}
#     intcoeffsplan::I  # getting coefficients of integrand to integrate
#     intcoeffs::Matrix{T}
#     intplan::P  # integration
# end

# function _compute!(s::ChebyshevSetup, Ω, D, phasepoints)
#     paduatransform!(s.phasecoeffs, s.paduaplan, phasepoints)
#     differentiate!(s.∂x, s.∂y, s.diffplan, s.phasecoeffs)
#     getintegrand!(s.intcoeffs, s.intcoeffsplan, Ω, phasepoints, s.phasecoeffs, s.∂x, s.∂y)
#     integrate(s.intcoeffs, s.intplan)
# end

end  # module ChebyshevImplementation
