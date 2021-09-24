module Chebyshev

using ApproxFunOrthogonalPolynomials
using FastTransforms: PaduaTransformPlan, plan_paduatransform!
using FastTransforms: paduavalsmat, trianglecfsvec!
using ApproxFunBase: TransformPlan, ITransformPlan, plan_transform, plan_itransform

using BlockBandedMatrices: BandedBlockBandedMatrix
using BlockArrays: BlockRange
using StaticArrays: SVector
using SparseArrays: sparse

using LinearAlgebra: mul!, rmul!, dot

getdegree(coeffnum) = (sqrt(1 + 8coeffnum) - 3) / 2
getcoeffnum(degree) = (degree + 1) * (degree + 2) ÷ 2

# Don't confuse with function getpointnum
# This function gets the number of points on the full 2D Chebyshev grid with polynomials,
# not just on the Padua grid
getfullpointnum(degree) = (degree + 1)^2

# happens to be same number since padua grid is halved
getpaduanum(degree) = getcoeffnum(degree)
nextpaduanum(N) = getpaduanum(ceil(Int, getdegree(N)))

function checkpaduanum(paduanum)
    @argcheck !isinteger(getdegree(paduanum)) "number of Padua points or coeffs must be a triangle number 1, 3, 6, 10, 15..."
end

struct PaduaSetup{T, PTP}
    plan::PTP
    coeffs::Matrix{T}
end

function PaduaSetup{T}(D, degree)
    N = getpaduanum(degree)

	# Val{true / false} indicates if its lexigraphical (i.e., x, y) or reverse (y, x)
	# Val{false} is the setting ApproxFun uses, which is why it's used here, too
    # Why is there a bang (!) at the end? It's not in-place.
	plan = plan_paduatransform!(T, paduanum, Val{false})
    coeffs = Matrix{T}(undef, N, D)

    PaduaSetup{T, typeof(plan)}(plan, coeffs)
end

function paduatransform!(
    out::AbstractVector{T},
    v::AbstractVector{T},
    P::PaduaTransformPlan
) where T
    axes(out, 1) == axes(v, 1) || error("axes of input and output vectors must match")
    N = length(v)
    checkpaduanum(N)
    n = Int(getdegree(N))
    
    vals = paduavalsmat(P, v)
    tensorcfs = P.dctplan * vals
    
    m, l = size(tensorcfs)
    rmul!(tensorcfs, T(2) / (n * (n + 1)))
    rmul!(view(tensorcfs,1,:), 0.5)
    rmul!(view(tensorcfs,:,1), 0.5)
    rmul!(view(tensorcfs,m,:), 0.5)
    rmul!(view(tensorcfs,:,l), 0.5)
    
    trianglecfsvec!(out, P, tensorcfs)

    out
end

function paduatransform!(setup::PaduaSetup, values)
    @argcheck axes(setup.coeffs) == axes(values) "axes of value array and preallocated coefficient array do not match"

    @views for i in axes(values, 2)
        paduatransform!(setup.coeffs[:, i], values[:, i], setup.plan)
    end

    return setup.coeffs
end

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

	D1raw = Derivative(CC, [1,0])[BlockRange(1:degree), BlockRange(1:degree+1)]
	UCtoUU = Conversion(UC, UU)[BlockRange(1:degree), BlockRange(1:degree)]
    Dx = sparse(UCtoUU * D1raw)

	D2raw = Derivative(CC, [0,1])[BlockRange(1:degree), BlockRange(1:degree+1)]
	CUtoUU = Conversion(CU, UU)[BlockRange(1:degree), BlockRange(1:degree)]
	Dy = sparse(CUtoUU * D2raw)

    # differentiating reduces degree by one
    coeffnum = getcoeffnum(degree - 1)

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

function ChebyshevSetup{T}(::AbstractMatrix, D::Integer, N::Integer)
    checkpaduanum(N)
    degree = getdegree(N)

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

"""
    getpaduapoints([::Type{T}, ]n::Integer)::Matrix{T} where {T <: Real}

returns padua points corresponding to degree `n` Chebyshev polynomial on square `0..1 × 0..1`
as matrix with element type `T`. Each row represens a point.
"""
function getpaduapoints(::Type{T}, n::Integer)::Matrix{T} where T <: Real
    paduanum = getpaduanum(n)
    out = Matrix{T}(undef, paduanum, 2)
    m = 0
    delta = 0
    NN = fld(n + 2, 2)
    @inbounds for k = n:-1:0
        if isodd(n)
            delta = mod(k, 2)
        end
        @inbounds for j = NN+delta:-1:1
            m += 1

            out[m, 1] = (sinpi(T(k) / T(n) - T(0.5)) + 1) / 2

            a = isodd(n - k) ? 1 : 2
            out[m, 2] = (sinpi(T(2j - a) / T(n + 1) - T(0.5)) + 1) / 2
        end
    end
    return out
end

getpaduapoints(n::Integer) = getpaduapoints(Float64, n)

getpoints(pinv::SecondPoincareInvariant{<:Any, T}) where T = getpaduapoints(T, ceil(Int, getdegree(pinv.N)))

end  # module Chebyshev
