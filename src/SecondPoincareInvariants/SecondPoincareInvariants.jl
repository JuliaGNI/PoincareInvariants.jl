@reexport module SecondPoincareInvariants

# Callable is Union{Function, Type}
using Base: Callable

using ..PoincareInvariants: AbstractPoincareInvariant, @argcheck
import ..PoincareInvariants: compute!, getpoints, getpointnum, getdim, getform

export SecondPoincareInvariant, PI2

## SecondPoincareInvariant ##

struct SecondPoincareInvariant{
    T,  # phase space and return type
    ΩT <: Union{Callable, AbstractMatrix},
    PS,
    P
} <: AbstractPoincareInvariant
    Ω::ΩT  # symplectic matrix or function returning one
    D::Int  # dimension of phase space
    pointspec::PS  # specifies how many points
    plan::P  # plan for chebyshev transform, differentiation, etc...
end

function SecondPoincareInvariant{T}(
    Ω::ΩT, D::Integer, N, P::Type=DEFAULT_PLAN_TYPE
) where {T, ΩT}
    ps = getpointspec(N, P)
    plan = P{T}(Ω, D, ps)
    SecondPoincareInvariant{T, ΩT, typeof(ps), typeof(plan)}(Ω, D, ps, plan)
end

function SecondPoincareInvariant{T}(Ω::ΩT, D::Integer, ps::PS, plan::P) where {T, ΩT, PS, P}
    SecondPoincareInvariant{T, ΩT, PS, P}(Ω, D, ps, plan)
end

const PI2 = SecondPoincareInvariant

getdim(pinv::SecondPoincareInvariant) = pinv.D
getform(pinv::SecondPoincareInvariant) = pinv.Ω

# Interface should define getpoints(f, T, N, P) method
getpoints(f::Function, pinv::SecondPoincareInvariant{T, Ω, NT, P}) where {T, Ω, NT, P} =
    getpoints(f, T, pinv.pointspec, P)
getpoints(pinv::SecondPoincareInvariant) = getpoints((x, y) -> (x, y), pinv)

getpointnum(pinv::SecondPoincareInvariant{T, Ω, PS, P}) where {T, Ω, PS, P} =
    getpointnum(pinv.pointspec, P)

getpointspec(N, P) = getpointnum(N, P)

## Implementation(s)

include("Chebyshev.jl")
include("FiniteDifferences.jl")

using .Chebyshev: ChebyshevPlan
using .FiniteDifferences: FiniteDiffPlan

const DEFAULT_PLAN_TYPE = ChebyshevPlan

end  # module SecondPoincareInvariants
