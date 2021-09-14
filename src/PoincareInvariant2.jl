get_uu_coeff_num(n) = n * (n + 1) ÷ 2

get_uu_point_num(n) = n^2

struct PoincareInvariant2{
	N,  # phase space dimension
	T <: Number,  # phase space type
	PTP <: PaduaTransformPlan,
	DBBB <: BandedBlockBandedMatrix,
	CBBB <: BandedBlockBandedMatrix,
	UUTP <: TransformPlan,
	UUITP <: ITransformPlan
} <: AbstractPoincareInvariant
	n::Int
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

function PoincareInvariant2{N, T}(padua_num::Integer) where {N, T}
	# n such that (n + 1) * (n + 2) ÷ 2 == padua_num
	# padua coefficients on upper triangular of matrix of size (n + 1) × (n + 1)
	check_padua_num(padua_num)
	n = Int(get_n(padua_num))

	# make plan for padua transform
	# Val{true / false} indicates if its lexigraphical (i.e., x, y) or reverse (y, x)
	# Val{false} is the setting ApproxFun uses, which is why it's used here, too
	padua_plan = plan_paduatransform!(T, padua_num, Val{false})

	# preallocate array for coefficients
	cc_coeffs = Matrix{T}(undef, padua_num, N)

	CC = Chebyshev(0..1)         ⊗ Chebyshev(0..1)
	UC = Ultraspherical(1, 0..1) ⊗ Chebyshev(0..1)
	CU = Chebyshev(0..1)         ⊗ Ultraspherical(1, 0..1)
	UU = Ultraspherical(1, 0..1) ⊗ Ultraspherical(1, 0..1)

	# preallocate Operators
	#
	# Each block represents a diagonal in the coefficient matrix
	# Differentiating removes highest order coefficients
	# hence we go from n + 1 blocks to n blocks
	#
	# differentiation converts to chebyshev basis of the second kind,
	# so we convert everything to chebyshev polynomials of the second kind
	D1 = Derivative(CC, [1,0])[BlockRange(1:n), BlockRange(1:n+1)]
	UCtoUU = Conversion(UC, UU)[BlockRange(1:n), BlockRange(1:n)]
	D1toUU = UCtoUU * D1

	D2 = Derivative(CC, [0,1])[BlockRange(1:n), BlockRange(1:n+1)]
	CUtoUU = Conversion(CU, UU)[BlockRange(1:n), BlockRange(1:n)]
	D2toUU = CUtoUU * D2

	# truncate highest order coefficients,
	# so number of coeffs matches number after differentiating
	CCtoUU = Conversion(CC, UU)[BlockRange(1:n), BlockRange(1:n+1)]
	
	uu_coeff_num = get_uu_coeff_num(n)
	uu_point_num = get_uu_point_num(n)
	
	uu_coeffs    = Matrix{T}(undef, uu_coeff_num, N)
	uu_d1_coeffs = Matrix{T}(undef, uu_coeff_num, N)
	uu_d2_coeffs = Matrix{T}(undef, uu_coeff_num, N)
	
	@assert size(D1toUU) == size(D2toUU) == size(CCtoUU) == (uu_coeff_num, padua_num)
	@assert typeof(D1toUU) == typeof(D2toUU) <: BandedBlockBandedMatrix

	uu_points = points(UU, uu_point_num) .|> SVector{2, T}
	uu_plan = plan_transform(UU, T, uu_point_num)  # values to coefficients
	uu_iplan = plan_itransform(UU, T, uu_coeff_num)  # coefficients to values
	
	uu_vals    = Matrix{T}(undef, uu_point_num, N)
	uu_d1_vals = Matrix{T}(undef, uu_point_num, N)
	uu_d2_vals = Matrix{T}(undef, uu_point_num, N)
		
	uu_I_vals = Vector{T}(undef, uu_point_num)

	UUIntegral = DefiniteIntegral(UU)[1:1, 1:uu_coeff_num]
		
	uu_I_coeffs = Vector{T}(undef, uu_coeff_num)

	PoincareInvariant2{N, T}(
		n, padua_plan, cc_coeffs, D1toUU, D2toUU, CCtoUU,
		uu_coeffs, uu_d1_coeffs, uu_d2_coeffs,
		uu_points, uu_iplan, uu_vals, uu_d1_vals, uu_d2_vals,
		uu_I_vals, uu_plan, uu_I_coeffs, UUIntegral
	)
end

function PoincareInvariant2{N, T}(
	n, padua_plan::PTP, cc_coeffs, D1toUU::DBBB, D2toUU::DBBB, CCtoUU::CBBB,
	uu_coeffs, uu_d1_coeffs, uu_d2_coeffs,
	uu_points, uu_iplan::UUITP, uu_vals, uu_d1_vals, uu_d2_vals,
	uu_I_vals, uu_plan::UUTP, uu_I_coeffs, UUIntegral
) where {N, T, PTP, DBBB, CBBB, UUITP, UUTP}
	PoincareInvariant2{N, T, PTP, DBBB, CBBB, UUTP, UUITP}(
		n, padua_plan, cc_coeffs, D1toUU, D2toUU, CCtoUU,
		uu_coeffs, uu_d1_coeffs, uu_d2_coeffs,
		uu_points, uu_iplan, uu_vals, uu_d1_vals, uu_d2_vals,
		uu_I_vals, uu_plan, uu_I_coeffs, UUIntegral
	)
end

function checkΩ(Ω::AbstractMatrix, N)
	size(Ω) == (N, N) || throw(ArgumentError(
		"Ω must be a $N × $N AbstractMatrix, not $(size(Ω, 1)) × $(size(Ω, 2))"
	))
end

function check_phase_points(phase_points::AbstractMatrix, n, N)
	num = get_padua_num(n)
	size(phase_points) == (num, N) || throw(ArgumentError(string(
		"phase_points must be a $num × $N AbstractMatrix,",
		"not $(size(phase_points, 1)) × $(size(phase_points, 2))"
	)))
end

function compute(
    pinv::PoincareInvariant2{N, T},
    phase_points::AbstractMatrix{T},
    Ω::AbstractMatrix
) where {N, T}
    checkΩ(Ω, N)
    check_phase_points(phase_points, pinv.n, N)

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

    for i in 1:get_uu_point_num(pinv.n)
        pinv.uu_I_vals[i] = dot(pinv.uu_d2_vals[i, :], Ω, pinv.uu_d1_vals[i, :])
    end

    # Convert to coefficients and calculate integral
    pinv.uu_I_coeffs .= pinv.uu_plan * copy(pinv.uu_I_vals)
    dot(pinv.UUIntegral, pinv.uu_I_coeffs)
end
