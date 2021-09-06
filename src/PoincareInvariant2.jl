struct PoincareInvariant2{
    N,  # phase space dimension
    T <: Number,  # phase space type
    PTP <: PaduaTransformPlan
} <: AbstractPoincareInvariant
    padua_plan::PTP
    # phase_coeffs::NTuple{N, Vector{T}}
end

function PoincareInvariant2(phase_points::NTuple{N, <:AbstractVector{T}}) where {N, T}
    #Assumes all have same length
    point_num = length(first(phase_points))

    # make plan for padua transform
    # Val{true / false} indicates if its lexigraphical (i.e., x, y) or reverse (y, x)
    padua_plan = plan_paduatransform!(T, point_num, Val{false})

    PoincareInvariant2{N, T}(padua_plan)
end

function PoincareInvariant2{N, T}(padua_plan) where {N, T}
    PoincareInvariant2{N, T, typeof(padua_plan)}(padua_plan)
end

function calculate!(
    pinv::PoincareInvariant2{N},
    Ω::AbstractMatrix
) where {N}

end