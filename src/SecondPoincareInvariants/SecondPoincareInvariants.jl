@reexport module SecondPoincareInvariants

# Callable is Union{Function, Type}
using Base: Callable

using ..PoincareInvariants: AbstractPoincareInvariant, @argcheck
import ..PoincareInvariants: compute!, getpoints, getpointnum, getdim

export SecondPoincareInvariant

## Implementations

include("ChebyshevImplementation.jl")

## SecondPoincareInvariant ##

# TODO: add isinplace check

struct SecondPoincareInvariant{
    T,  # phase space and return type
    ΩT <: Union{Callable, AbstractMatrix},
    P
} <: AbstractPoincareInvariant
    Ω::ΩT  # symplectic matrix or function returning one
    D::Int
    plan::P  # plan for chebyshev transform, differentiation, etc...
end

function SecondPoincareInvariant{T}(Ω::ΩT, D::Integer, N::Integer) where {T, ΩT <: AbstractMatrix}
    @argcheck size(Ω) == (D, D) "Ω must be a $D × $D matrix"
    plan = ChebyshevImplementation.ChebyshevPlan{T}(Ω, D, N)
    SecondPoincareInvariant{T, ΩT, typeof(plan)}(Ω, D, plan)
end

function SecondPoincareInvariant{T}(
    Ω::ΩT, D::Integer, N::Integer, ::Val{inplace}
) where {T, ΩT <: Callable, inplace}
    plan = ChebyshevImplementation.ChebyshevPlan{T}(Ω, D, N, Val(inplace))
    SecondPoincareInvariant{T, ΩT, typeof(plan)}(Ω, D, plan)
end

# Unexported convenience alias
const PI2 = SecondPoincareInvariant

getdim(pinv::SecondPoincareInvariant) = pinv.D

## Internal interface to implementation ##
compute!(pinv::SecondPoincareInvariant, phasepoints::AbstractMatrix, t, p) =
    compute!(pinv.plan, pinv.Ω, phasepoints, t, p)
compute!(pinv::SecondPoincareInvariant, phasepoints::AbstractMatrix) =
    compute!(pinv.plan, pinv.Ω, phasepoints)

getpoints(pinv::SecondPoincareInvariant) = getpoints(pinv.plan)
getpoints(f::Function, pinv::SecondPoincareInvariant) = getpoints(f, pinv.plan)
getpointnum(pinv::SecondPoincareInvariant) = getpointnum(pinv.plan)

end  # module SecondPoincareInvariants
