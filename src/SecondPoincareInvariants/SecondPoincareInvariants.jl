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
    NT,
    P
} <: AbstractPoincareInvariant
    Ω::ΩT  # symplectic matrix or function returning one
    D::Int  # dimension of phase space
    N::NT  # specifies how many points
    plan::P  # plan for chebyshev transform, differentiation, etc...
end

function SecondPoincareInvariant{T}(
    Ω::ΩT, D::Integer, N::NT, P::Type=DEFAULT_PLAN_TYPE
) where {T, ΩT, NT}
    plan = P{T}(Ω, D, N)
    SecondPoincareInvariant{T, ΩT, NT, typeof(plan)}(Ω, D, N, plan)
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

# Interface should define getpoints(f, T, N, P) method
getpoints(f::Function, pinv::SecondPoincareInvariant{T, Ω, NT, P}) where {T, Ω, NT, P} =
    getpoints(f, T, pinv.N, P)
getpoints(pinv::SecondPoincareInvariant) = getpoints((x, y) -> (x, y), pinv)

getpointnum(pinv::SecondPoincareInvariant{T, Ω, P}) where {T, Ω, NT, P} =
    getpointnum(pinv.N, P)

## Implementation(s)

include("Chebyshev.jl")
include("Trapezoidal.jl")

const DEFAULT_PLAN_TYPE = Chebyshev.ChebyshevPlan

end  # module SecondPoincareInvariants
