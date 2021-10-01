"""
    ChebyshevImplementation

implementation of computation of second Poincare invariant by approximating surface with
Chebyshev polynomials
"""
module ChebyshevImplementation

using ...PoincareInvariants: @argcheck

using ApproxFunOrthogonalPolynomials
using ApproxFunBase: TransformPlan, ITransformPlan, plan_transform, plan_itransform

using BlockArrays: BlockRange
using StaticArrays: SVector
using SparseArrays: sparse

using LinearAlgebra: mul!, rmul!, dot

include("Padua.jl")
import .Padua: paduatransform!

using .Padua


# This function gets the number of points on the full 2D Chebyshev grid with polynomials,
# not just on the Padua grid
getfullpointnum(degree) = (degree + 1)^2

## Padua Transform ##

struct PaduaSetup{T, PTP}
    plan::PTP
    coeffs::Matrix{T}
end

function PaduaSetup{T}(D, degree) where T
    N = getpaduanum(degree)

	plan = PaduaTransformPlan{T}(degree)
    coeffs = Matrix{T}(undef, N, D)

    PaduaSetup{T, typeof(plan)}(plan, coeffs)
end

function paduatransform!(setup::PaduaSetup, values, args...)
    paduatransform!(setup.coeffs, setup.plan, values, args...)
end

## Differentiation ##

const CC = Chebyshev(-1..1)         ⊗ Chebyshev(-1..1)
const UC = Ultraspherical(1, -1..1) ⊗ Chebyshev(-1..1)
const CU = Chebyshev(-1..1)         ⊗ Ultraspherical(1, -1..1)
const UU = Ultraspherical(1, -1..1) ⊗ Ultraspherical(1, -1..1)

struct DiffSetup{T, DM}
    Dx::DM
    Dy::DM
    ∂xcoeffs::Matrix{T}
    ∂ycoeffs::Matrix{T}
end

function DiffSetup{T}(D, degree) where T
	N = getpaduanum(degree)

	# Each block represents a diagonal in the coefficient matrix
	# Differentiating removes highest order coefficients
	# hence we go from n + 1 blocks to n blocks
	#
	# differentiation converts to chebyshev basis of the second kind,
	# so we convert other axis to chebyshev of second kind, too

	D1UC = Derivative(CC, [1,0])[BlockRange(1:degree), BlockRange(1:degree+1)]
	UCtoUU = Conversion(UC, UU)[BlockRange(1:degree), BlockRange(1:degree)]
    Dx = sparse(UCtoUU * D1UC)

	D2CU = Derivative(CC, [0,1])[BlockRange(1:degree), BlockRange(1:degree+1)]
	CUtoUU = Conversion(CU, UU)[BlockRange(1:degree), BlockRange(1:degree)]
	Dy = sparse(CUtoUU * D2CU)

    # differentiating reduces degree by one
    coeffnum = getpaduanum(degree - 1)

    ∂xcoeffs = Matrix{T}(undef, coeffnum, D)
	∂ycoeffs = Matrix{T}(undef, coeffnum, D)

    typeof(Dx) == typeof(Dy) || error()

    DiffSetup{T, typeof(Dx)}(Dx, Dy, ∂xcoeffs, ∂ycoeffs)
end

function differentiate!(setup::DiffSetup, coeffs::AbstractMatrix)
    mul!(setup.∂xcoeffs, setup.Dx, coeffs)
    mul!(setup.∂ycoeffs, setup.Dy, coeffs)
    return setup.∂xcoeffs, setup.∂ycoeffs
end

struct C12Setup{T, CM}
    C12::CM
    coeffs::Matrix{T}
end

function C12Setup{T}(D, degree) where T
    N = getpaduanum(degree)

    C12 = sparse(Conversion(CC, UU)[BlockRange(1:degree), BlockRange(1:degree+1)])

    # truncate degree so degree of derivative is the same as converted polynomial
    coeffnum = getcoeffnum(degree - 1)
    
    coeffs = Matrix{T}(undef, coeffnum, D)

    C12Setup{T, typeof(C12)}(C12, coeffs)
end

# convert from chebyshev polynomials of first to second kind
C12convert!(setup::C12Setup, coeffs) = mul!(setup.coeffs, setup.C12, coeffs)

struct ConstΩIntSetup{T, C12T <: C12Setup{T}, ITP, TP}
    C12::C12T
    iplan::ITP
    phasevals::Matrix{T}
    ∂xvals::Matrix{T}
    ∂yvals::Matrix{T}
    integrandvals::Vector{T}
    plan::TP
    integrandcoeffs::Vector{T}
    Integral::Matrix{T}
end

function ConstΩIntSetup{T}(D, degree) where T
    coeffnum = getcoeffnum(degree - 1)
    pointnum = getfullpointnum(degree - 1)

    C12 = C12Setup{T}(D, degree)
    iplan = plan_itransform(UU, T, coeffnum)

    phasevals = Matrix{T}(undef, pointnum, D)
    ∂xvals = Matrix{T}(undef, pointnum, D)
    ∂yvals = Matrix{T}(undef, pointnum, D)

    integrandvals = Vector{T}(undef, pointnum)

    plan = plan_transform(UU, T, pointnum)

    integrandcoeffs = Vector{T}(undef, coeffnum)

    Integral = DefiniteIntegral(UU)[1:1, 1:coeffnum]

    ConstΩIntSetup{T, typeof(C12), typeof(iplan), typeof(plan)}(
        C12, iplan, phasevals, ∂xvals, ∂yvals, integrandvals, plan, integrandcoeffs, Integral
    )
end

function integrate!(
    setup::ConstΩIntSetup, Ω::AbstractMatrix, D::Integer,
    phasecoeffs, ∂xcoeffs, ∂ycoeffs, t, p
)
    coeffs = C12convert!(setup.C12, phasecoeffs)
    
    # Convert back to values from coefficients
    @views for i in 1:D
        copyto!(setup.phasevals[:, i], coeffs[:, i])
        setup.iplan * setup.phasevals[:, i]

        copyto!(setup.∂xvals[:, i], ∂xcoeffs[:, i])
        setup.iplan * setup.∂xvals[:, i]

        copyto!(setup.∂yvals[:, i], ∂ycoeffs[:, i])
        setup.iplan * setup.∂yvals[:, i]
    end

    # Calculate integrand on full chebyshev grid
    @views for i in 1:length(setup.integrandvals)
        setup.integrandvals[i] = dot(setup.∂yvals[i, :], Ω, setup.∂xvals[i, :])
    end

    # Convert to coefficients and calculate integral
    copyto!(setup.integrandcoeffs, setup.integrandvals)
    setup.iplan * setup.integrandcoeffs
    dot(setup.Integral, setup.integrandcoeffs)
end

struct ChebyshevSetup{T, I}
    degree::Int
    padua::PaduaSetup
    diff::DiffSetup
    int::I
end

function ChebyshevSetup{T}(::AbstractMatrix, D::Integer, N::Integer) where T
    degree = ceil(Int, getdegree(N))

    padua = PaduaSetup{T}(D, degree)
    diff = DiffSetup{T}(D, degree)
    int = ConstΩIntSetup{T}(D, degree)

    ChebyshevSetup{T, typeof(int)}(degree, padua, diff, int)
end

function _compute(setup::ChebyshevSetup, Ω, D, phasepoints, t, p)
	# Convert to basis of Chebyshev polynomials of the first kind
	coeffs = paduatransform!(setup.padua, phasepoints)
	
	# Get coefficients of derivatives in basis of Chebyshev polynomials of the second kind
	∂x, ∂y = differentiate!(setup.diff, coeffs)

	# Calculate ∑ ∫ Ωᵢⱼ(z(x, y), p, t)⋅(∂z/∂x)ᵢ⋅(∂z/∂y)ⱼ dx dy
	integrate!(setup.int, Ω, D, C1coeffs, C2∂x, C2∂y, t, p)
end

getpoints(setup::ChebyshevSetup) where T = getpaduapoints(T, setup.degree)

end  # module ChebyshevImplementation
