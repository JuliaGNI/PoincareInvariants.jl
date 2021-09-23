struct Setup2{
	T <: Number,  # phase space type
	PTP <: PaduaTransformPlan,
	DM <: AbstractMatrix,
	CM <: AbstractMatrix,
	UUTP <: TransformPlan,
	UUITP <: ITransformPlan
} <: AbstractPoincareInvariant
	paduaplan::PTP
	C1coeffs::AbstractMatrix{T}
	DutoC2::DM
	DvtoC2::DM
	C1toC2::CM
	C2coeffs::Matrix{T}
	C2Ducoeffs::Matrix{T}
	C2Dvcoeffs::Matrix{T}
	C2iplan::UUITP
	C2vals::Matrix{T}
	C2Duvals::Matrix{T}
	C2Dvvals::Matrix{T}
	C2Ivals::Vector{T}
	C2plan::UUTP
	C2Icoeffs::Vector{T}
	C2Integral::Matrix{T}
end

function Setup2(D, ::Type{T}, ::AbstractMatrix, paduanum) where T
	checkpaduanum(paduanum)
	degree = Int(getdegree(paduanum))

	# make plan for Padua transform
	# Val{true / false} indicates if its lexigraphical (i.e., x, y) or reverse (y, x)
	# Val{false} is the setting ApproxFun uses, which is why it's used here, too
    # Why is there a bang (!) at the end? It's not in-place.
	paduaplan = plan_paduatransform!(T, paduanum, Val{false})

	# preallocate array for coefficients
	C1coeffs = Matrix{T}(undef, paduanum, D)

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
	Du = Derivative(CC, [1,0])[BlockRange(1:degree), BlockRange(1:degree+1)]
	UCtoC2 = Conversion(UC, UU)[BlockRange(1:degree), BlockRange(1:degree)]
	DutoC2 = sparse(UCtoC2 * Du)

	Dv = Derivative(CC, [0,1])[BlockRange(1:degree), BlockRange(1:degree+1)]
	CUtoC2 = Conversion(CU, UU)[BlockRange(1:degree), BlockRange(1:degree)]
	DvtoC2 = sparse(CUtoC2 * Dv)

	# truncate highest order coefficients,
	# so number of coeffs matches number after differentiating
	C1toC2 = sparse(Conversion(CC, UU)[BlockRange(1:degree), BlockRange(1:degree+1)])
	
	C2coeffnum = getcoeffnum(degree - 1)
	C2pointnum = getpointnum(degree - 1)
	
	C2coeffs   = Matrix{T}(undef, C2coeffnum, D)
	C2Ducoeffs = Matrix{T}(undef, C2coeffnum, D)
	C2Dvcoeffs = Matrix{T}(undef, C2coeffnum, D)

	C2plan = plan_transform(UU, T, C2pointnum)  # values to coefficients
	C2iplan = plan_itransform(UU, T, C2coeffnum)  # coefficients to values

	C2vals   = Matrix{T}(undef, C2pointnum, D)
	C2Duvals = Matrix{T}(undef, C2pointnum, D)
	C2Dvvals = Matrix{T}(undef, C2pointnum, D)

	C2Ivals = Vector{T}(undef, C2pointnum)

	C2Integral = DefiniteIntegral(UU)[1:1, 1:C2coeffnum]

	C2Icoeffs = Vector{T}(undef, C2coeffnum)

	Setup2{
		T, typeof(paduaplan), typeof(DutoC2), typeof(C1toC2),
		typeof(C2plan), typeof(C2iplan)
	}(
		paduaplan, C1coeffs, DutoC2, DvtoC2, C1toC2,
		C2coeffs, C2Ducoeffs, C2Dvcoeffs,
		C2iplan, C2vals, C2Duvals, C2Dvvals,
		C2Ivals, C2plan, C2Icoeffs, C2Integral
	)
end

function checkphasepoints(phasepoints::AbstractMatrix, N, D)
	size(phasepoints) == (N, D) || throw(ArgumentError(string(
		"phasepoints must be a $N × $D and not a",
		"$(size(phasepoints, 1)) × $(size(phasepoints, 2)) AbstractMatrix"
	)))
end

function compute(
    pinv::SecondPoincareInvariant{D, T, ΩT, S},
    phasepoints::AbstractMatrix
) where {D, T, ΩT, S <: Setup2}
	checkphasepoints(phasepoints, pinv.N, D)

	setup = pinv.setup

    # Transform to Chebyshev basis via Padua transform
    for i in 1:D
        setup.C1coeffs[:, i] .= setup.paduaplan * copy(phasepoints[:, i])
    end

    # Calculate derivatives and convert to chebyshev basis of second kind
    mul!(  setup.C2coeffs, setup.C1toC2, setup.C1coeffs)
    mul!(setup.C2Ducoeffs, setup.DutoC2, setup.C1coeffs)
    mul!(setup.C2Dvcoeffs, setup.DvtoC2, setup.C1coeffs)

    # Convert back to values to get dot product with Ω
    for i in 1:D
        setup.C2vals[:, i]   .= setup.C2iplan * copy(setup.C2coeffs[:, i])
        setup.C2Duvals[:, i] .= setup.C2iplan * copy(setup.C2Ducoeffs[:, i])
        setup.C2Dvvals[:, i] .= setup.C2iplan * copy(setup.C2Dvcoeffs[:, i])
    end

    for i in 1:Int(getpointnum(getdegree(pinv.N) - 1))
        setup.C2Ivals[i] = dot(setup.C2Dvvals[i, :], pinv.Ω, setup.C2Duvals[i, :])
    end

    # Convert to coefficients and calculate integral
    setup.C2Icoeffs .= setup.C2plan * copy(setup.C2Ivals)
    dot(setup.C2Integral, setup.C2Icoeffs)
end
