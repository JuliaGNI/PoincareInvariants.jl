module PoincareInvariants

# For the Padua transforms
using FastTransforms: PaduaTransformPlan, plan_paduatransform!

# to represent points
using StaticArrays

# Callable is Union{Function, Type}
using Base: Callable

export AbstractPoincareInvariant, PoincareInvariant2

const PARAM_TYPE = Float64

# Do we want PoincareInvariant2 <: AbstractPoincareInvariant2 <: AbstractPoincareInvariant ?
# I guess AbstractPoincareInvariant2 would be worth it if we had different types of
# PoincareInvariant2, like CanonicalPoincareInvariant2 may be
abstract type AbstractPoincareInvariant end

# get minimum allowed padua point number ≥ mn
function get_min_padua_num(mn)
    n = Int(cld(-3 + sqrt(1 + 8mn), 2))
    return (n + 1) * (n + 2) ÷ 2
end

# Need to wrap paduapoints function since it returns a matrix of the points
# and a Vector of SVector's is easier to work with
function _get_padua_points(::Type{T}, point_num::Integer)::Vector{SVector{2, T}} where T
    mat = paduapoints(T, Int(cld(-3 + sqrt(1 + 8point_num), 2)))
    map(eachrow(mat)) do row
        # rescale from -1..1 to 0..1
        (SVector{2}(row) .+ 1) ./ 2
    end
end

# Should the user be free to choose parameter space type?
struct PoincareInvariant2{
    N,  # phase space dimension
    T <: Number,  # phase space type
    PF <: Callable,
    PTP <: PaduaTransformPlan
} <: AbstractPoincareInvariant
    param_func::PF
    point_num::Int  # make this part of type?
    padua_points::Vector{SVector{2, PARAM_TYPE}}
    phase_points::NTuple{N, Vector{T}}
    padua_plan::PTP
end

# This is the constructor the user should use
function PoincareInvariant2{N}(param_func, min_point_num) where N
    # only certain number of padua points possible
    point_num = get_min_padua_num(min_point_num)

    # allocate phase_points
    T = SVector{2, PARAM_TYPE}(0, 0) |> param_func |> eltype
    phase_points = ntuple(_ -> Vector{T}(undef, point_num), Val(N))

    # get padua points
    padua_points = _get_padua_points(PARAM_TYPE, point_num)

    # make plan for padua transform
    # Val{true / false} indicates if its lexigraphical (i.e., x, y) or reverse (y, x)
    padua_plan = plan_paduatransform!(T, point_num, Val{false})

    PoincareInvariant2{N, T}(param_func, point_num, padua_points, phase_points, padua_plan)
end

# fill in types
function PoincareInvariant2{N, T}(
    param_func, point_num, padua_points, phase_points, padua_plan
) where {N, T}
    PoincareInvariant2{N, T, typeof(param_func), typeof(padua_plan)}(
        param_func, point_num, padua_points, phase_points, padua_plan
    )
end

# calculates phase points using param_func and outputs a vector of values for each
# dimension of phase space (as opposed to vector of SVectors)
#
# Type annotations added to force specialisation on function
#
# TODO: add @generated version with @nexprs
# TODO: Should generated version also use known point number?
function get_phase_points!(param_func::F, out, param_points, ::Val{N}) where {F, N}
    N::Int  # base Julia ntuple.jl does this, too
    
    for (i, pa_pnt) in enumerate(param_points)
        ph_pnt = param_func(pa_pnt)
        for j in 1:N
            out[j][i] = ph_pnt[j]
        end
    end

    out
end

end  # module
