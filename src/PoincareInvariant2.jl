struct PoincareInvariant2{
    N,  # phase space dimension
    T <: Number,  # phase space type
    PTP <: PaduaTransformPlan
} <: AbstractPoincareInvariant
    padua_plan::PTP
    # param_coeffs::NTuple{N, Vector{T}}
end

function PoincareInvariant2(phase_points::AbstractMatrix{T}) where {T}
    #Assumes all have same length
    point_num, N = size(phase_points)

    # make plan for padua transform
    # Val{true / false} indicates if its lexigraphical (i.e., x, y) or reverse (y, x)
    padua_plan = plan_paduatransform!(T, point_num, Val{false})

    PoincareInvariant2{N, T}(padua_plan)
end

function PoincareInvariant2{N, T}(padua_plan) where {N, T}
    PoincareInvariant2{N, T, typeof(padua_plan)}(padua_plan)
end

const CC = Chebyshev(0..1)         ⊗ Chebyshev(0..1)
const UC = Ultraspherical(1, 0..1) ⊗ Chebyshev(0..1)
const CU = Chebyshev(0..1)         ⊗ Ultraspherical(1, 0..1)
const UU = Ultraspherical(1, 0..1) ⊗ Ultraspherical(1, 0..1)

const CCtoUU = Conversion(CC, UU)

const D1 = Derivative(CC, [1,0])
const UCtoUU = Conversion(UC, UU)

const D2 = Derivative(CC, [0,1])
const CUtoUU = Conversion(CU, UU)

const UUDefInt = DefiniteIntegral(UU)

function compute_uu_phase_vals_and_derivs(phase_vals, plan)
    np = length(phase_vals)
    
    # phase vals to coefficients in C^2
    param_coeffs = similar(phase_vals)
    padua_transform!(param_coeffs, phase_vals, plan)
    pf = Fun(CC, param_coeffs)

    # differentiate and convert to U^2
    Dpf1 = UCtoUU * (D1 * pf)
    Dpf2 = CUtoUU * (D2 * pf)

    # itransform to vals on U^2 grid
    # derivative has one fewer point. So, this drops highest order coefficient
    uu_param_coeffs = resize!(collect((CCtoUU * pf).coefficients), ncoefficients(Dpf1))
    uu_phase_vals = itransform(UU, uu_param_coeffs)

    return uu_phase_vals, values(Dpf1), values(Dpf2)
end

compute_integral(Is) = Number(UUDefInt * Fun(UU, transform(UU, Is)))

function compute(
    pinv::PoincareInvariant2{N},
    phase_points::AbstractMatrix{T},
    Ω::Callable
) where {N, T}
    np = size(phase_points, 1)

    vals, d1, d2  = compute_uu_phase_vals_and_derivs(phase_points[:, 1], pinv.padua_plan)

    uu_phase_points = Matrix{T}(undef, length(vals), N)
    uu_d1s = Matrix{T}(undef, length(d1), N)
    uu_d2s = Matrix{T}(undef, length(d2), N)

    uu_phase_points[:, 1] .= vals
    uu_d1s[:, 1] .= d1
    uu_d2s[:, 1] .= d2

    for i in 2:N
        vals, d1, d2  = compute_uu_phase_vals_and_derivs(phase_points[:, i], pinv.padua_plan)
        
        uu_phase_points[:, i] .= vals
        uu_d1s[:, i] .= d1
        uu_d2s[:, i] .= d2 
    end

    uu_Is = similar(uu_phase_points, size(uu_phase_points, 1))

    # compute integrands of integral invariants at all points
    for (i, point) in enumerate(eachrow(uu_phase_points))
        uu_Is[i] = dot(uu_d2s[i, :], Ω(point), uu_d1s[i, :])
    end

    # compute integral invariant
    return compute_integral(uu_Is)
end
