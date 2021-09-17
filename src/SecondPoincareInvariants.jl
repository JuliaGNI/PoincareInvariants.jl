using ApproxFunOrthogonalPolynomials
using FastTransforms: PaduaTransformPlan, plan_paduatransform!
using ApproxFunBase: TransformPlan, ITransformPlan, plan_transform, plan_itransform

using BlockBandedMatrices: BandedBlockBandedMatrix
using BlockArrays: BlockRange
using StaticArrays: SVector

using LinearAlgebra: mul!, dot

# Callable is Union{Function, Type}
using Base: Callable

using ..PoincareInvariants: AbstractPoincareInvariant
import ..PoincareInvariants: compute

export SecondPoincareInvariant
export get_padua_points, next_padua_num

get_n(padua_num) = (sqrt(1 + 8padua_num) - 3) / 2

get_padua_num(degree) = (degree + 1) * (degree + 2) ÷ 2

next_padua_num(N) = get_padua_num(ceil(Int, get_n(N)))

function check_padua_num(padua_num)
    if !isinteger(get_n(padua_num))
        throw(ArgumentError("number of padua points must equal (n + 1) * (n + 2) ÷ 2"))
    end
end

include("padua.jl")

get_uu_coeff_num(degree) = degree * (degree + 1) ÷ 2

get_uu_point_num(degree) = degree^2

## SecondPoincareInvariant ##

struct SecondPoincareInvariant{
	N,  # phase space dimension
	T <: Number,  # phase space type
	PTP <: PaduaTransformPlan,
	DBBB <: BandedBlockBandedMatrix,
	CBBB <: BandedBlockBandedMatrix,
	UUTP <: TransformPlan,
	UUITP <: ITransformPlan
} <: AbstractPoincareInvariant
	degree::Int
	padua_plan::PTP
	cc_coeffs::AbstractMatrix{T}
	D1toUU::DBBB
	D2toUU::DBBB
	CCtoUU::CBBB
	uu_coeffs::Matrix{T}
	uu_d1_coeffs::Matrix{T}
	uu_d2_coeffs::Matrix{T}
	uu_points::Vector{SVector{2, T}}
	uu_iplan::UUITP
	uu_vals::Matrix{T}
	uu_d1_vals::Matrix{T}
	uu_d2_vals::Matrix{T}
	uu_I_vals::Vector{T}
	uu_plan::UUTP
	uu_I_coeffs::Vector{T}
	UUIntegral::Matrix{T}
end

function SecondPoincareInvariant{N, T}(padua_num::Integer) where {N, T}
	# n such that (n + 1) * (n + 2) ÷ 2 == padua_num
	# padua coefficients on upper triangular of matrix of size (n + 1) × (n + 1)
	check_padua_num(padua_num)
	degree = Int(get_n(padua_num))

	# make plan for padua transform
	# Val{true / false} indicates if its lexigraphical (i.e., x, y) or reverse (y, x)
	# Val{false} is the setting ApproxFun uses, which is why it's used here, too
    # Why is there a bang (!) at the end? It's not in-place.
	padua_plan = plan_paduatransform!(T, padua_num, Val{false})

	# preallocate array for coefficients
	cc_coeffs = Matrix{T}(undef, padua_num, N)

	CC = Chebyshev(-1..1)         ⊗ Chebyshev(-1..1)
	UC = Ultraspherical(1, -1..1) ⊗ Chebyshev(-1..1)
	CU = Chebyshev(-1..1)         ⊗ Ultraspherical(1, -1..1)
	UU = Ultraspherical(1, -1..1) ⊗ Ultraspherical(1, -1..1)

	# preallocate Operators
	#
	# Each block represents a diagonal in the coefficient matrix
	# Differentiating removes highest order coefficients
	# hence we go from n + 1 blocks to n blocks
	#
	# differentiation converts to chebyshev basis of the second kind,
	# so we convert everything to chebyshev polynomials of the second kind
	D1 = Derivative(CC, [1,0])[BlockRange(1:degree), BlockRange(1:degree+1)]
	UCtoUU = Conversion(UC, UU)[BlockRange(1:degree), BlockRange(1:degree)]
	D1toUU = UCtoUU * D1

	D2 = Derivative(CC, [0,1])[BlockRange(1:degree), BlockRange(1:degree+1)]
	CUtoUU = Conversion(CU, UU)[BlockRange(1:degree), BlockRange(1:degree)]
	D2toUU = CUtoUU * D2

	# truncate highest order coefficients,
	# so number of coeffs matches number after differentiating
	CCtoUU = Conversion(CC, UU)[BlockRange(1:degree), BlockRange(1:degree+1)]
	
	uu_coeff_num = get_uu_coeff_num(degree)
	uu_point_num = get_uu_point_num(degree)
	
	uu_coeffs    = Matrix{T}(undef, uu_coeff_num, N)
	uu_d1_coeffs = Matrix{T}(undef, uu_coeff_num, N)
	uu_d2_coeffs = Matrix{T}(undef, uu_coeff_num, N)

	uu_points = points(UU, uu_point_num) .|> SVector{2, T}
	uu_plan = plan_transform(UU, T, uu_point_num)  # values to coefficients
	uu_iplan = plan_itransform(UU, T, uu_coeff_num)  # coefficients to values

	uu_vals    = Matrix{T}(undef, uu_point_num, N)
	uu_d1_vals = Matrix{T}(undef, uu_point_num, N)
	uu_d2_vals = Matrix{T}(undef, uu_point_num, N)

	uu_I_vals = Vector{T}(undef, uu_point_num)

	UUIntegral = DefiniteIntegral(UU)[1:1, 1:uu_coeff_num]

	uu_I_coeffs = Vector{T}(undef, uu_coeff_num)

	SecondPoincareInvariant{N, T}(
		degree, padua_plan, cc_coeffs, D1toUU, D2toUU, CCtoUU,
		uu_coeffs, uu_d1_coeffs, uu_d2_coeffs,
		uu_points, uu_iplan, uu_vals, uu_d1_vals, uu_d2_vals,
		uu_I_vals, uu_plan, uu_I_coeffs, UUIntegral
	)
end

function SecondPoincareInvariant{N, T}(
	degree, padua_plan::PTP, cc_coeffs, D1toUU::DBBB, D2toUU::DBBB, CCtoUU::CBBB,
	uu_coeffs, uu_d1_coeffs, uu_d2_coeffs,
	uu_points, uu_iplan::UUITP, uu_vals, uu_d1_vals, uu_d2_vals,
	uu_I_vals, uu_plan::UUTP, uu_I_coeffs, UUIntegral
) where {N, T, PTP, DBBB, CBBB, UUITP, UUTP}
	SecondPoincareInvariant{N, T, PTP, DBBB, CBBB, UUTP, UUITP}(
		degree, padua_plan, cc_coeffs, D1toUU, D2toUU, CCtoUU,
		uu_coeffs, uu_d1_coeffs, uu_d2_coeffs,
		uu_points, uu_iplan, uu_vals, uu_d1_vals, uu_d2_vals,
		uu_I_vals, uu_plan, uu_I_coeffs, UUIntegral
	)
end

const PI2 = SecondPoincareInvariant

## compute ##

function checkΩ(Ω::AbstractMatrix, N)
	size(Ω) == (N, N) || throw(ArgumentError(
		"Ω must be a $N × $N AbstractMatrix, not $(size(Ω, 1)) × $(size(Ω, 2))"
	))
end

function check_phase_points(phase_points::AbstractMatrix, degree, N)
	num = get_padua_num(degree)
	size(phase_points) == (num, N) || throw(ArgumentError(string(
		"phase_points must be a $num × $N AbstractMatrix,",
		"not $(size(phase_points, 1)) × $(size(phase_points, 2))"
	)))
end

function compute(
    pinv::SecondPoincareInvariant{N, T},
    phase_points::AbstractMatrix,
    Ω::AbstractMatrix
) where {N, T}
    checkΩ(Ω, N)
    check_phase_points(phase_points, pinv.degree, N)

    # Transform to Chebyshev basis via Padua transform
    for i in 1:N
        pinv.cc_coeffs[:, i] .= pinv.padua_plan * copy(phase_points[:, i])
    end

    # Calculate derivatives and convert to chebyshev basis of second kind
    mul!(   pinv.uu_coeffs, pinv.CCtoUU, pinv.cc_coeffs)
    mul!(pinv.uu_d1_coeffs, pinv.D1toUU, pinv.cc_coeffs)
    mul!(pinv.uu_d2_coeffs, pinv.D2toUU, pinv.cc_coeffs)

    # Convert back to values to get dot product with Ω
    for i in 1:N
        pinv.uu_vals[:, i]    .= pinv.uu_iplan * copy(pinv.uu_coeffs[:, i])
        pinv.uu_d1_vals[:, i] .= pinv.uu_iplan * copy(pinv.uu_d1_coeffs[:, i])
        pinv.uu_d2_vals[:, i] .= pinv.uu_iplan * copy(pinv.uu_d2_coeffs[:, i])
    end

    for i in 1:get_uu_point_num(pinv.degree)
        pinv.uu_I_vals[i] = dot(pinv.uu_d2_vals[i, :], Ω, pinv.uu_d1_vals[i, :])
    end

    # Convert to coefficients and calculate integral
    pinv.uu_I_coeffs .= pinv.uu_plan * copy(pinv.uu_I_vals)
    dot(pinv.UUIntegral, pinv.uu_I_coeffs)
end
