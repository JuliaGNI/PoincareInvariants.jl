@reexport module SecondPoincareInvariants

# Callable is Union{Function, Type}
using Base: Callable

using ..PoincareInvariants: AbstractPoincareInvariant, @argcheck
import ..PoincareInvariants: compute!, getpoints, getpointnum, getdim, getform

export SecondPoincareInvariant

## SecondPoincareInvariant ##

struct SecondPoincareInvariant{
    T,  # phase space and return type
    ΩT <: Union{Callable, AbstractMatrix},
    P
} <: AbstractPoincareInvariant
    Ω::ΩT  # symplectic matrix or function returning one
    D::Int
    plan::P  # plan for chebyshev transform, differentiation, etc...
end

function SecondPoincareInvariant{T}(
    Ω::ΩT, D::Integer, N::Integer, P::Type=DEFAULT_PLAN_TYPE
) where {T, ΩT}
    plan = P{T}(Ω, D, N)
    SecondPoincareInvariant{T, ΩT, typeof(plan)}(Ω, D, plan)
end

# Unexported convenience alias
# TODO: move or export these aliases to top level module?
const PI2 = SecondPoincareInvariant

getdim(pinv::SecondPoincareInvariant) = pinv.D
getform(pinv::SecondPoincareInvariant) = pinv.Ω

## Internal interface to implementation ##

compute!(pinv::SecondPoincareInvariant, phasepoints, t, p) =
    compute!(pinv.plan, pinv.Ω, phasepoints, t, p)
compute!(pinv::SecondPoincareInvariant, phasepoints) =
    compute!(pinv.plan, pinv.Ω, phasepoints)

getpoints(pinv::SecondPoincareInvariant) = getpoints(pinv.plan)
getpoints(f::Function, pinv::SecondPoincareInvariant) = getpoints(f, pinv.plan)
getpointnum(pinv::SecondPoincareInvariant) = getpointnum(pinv.plan)

## Implementation(s)

include("Chebyshev.jl")
include("FiniteDifferences.jl")

const DEFAULT_PLAN_TYPE = Chebyshev.ChebyshevPlan

end  # module SecondPoincareInvariants
