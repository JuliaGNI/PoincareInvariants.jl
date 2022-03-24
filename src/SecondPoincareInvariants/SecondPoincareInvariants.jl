@reexport module SecondPoincareInvariants

import ..PoincareInvariants: compute!, getpoints, getpointspec, getdim, getform

using ..PoincareInvariants: AbstractPoincareInvariant, @argcheck

# Callable is Union{Function, Type}
using Base: Callable

export SecondPoincareInvariant, PI2

## SecondPoincareInvariant ##

struct SecondPoincareInvariant{
    T,  # phase space and return type
    D,  # dimension of phase space
    ΩT <: Union{Callable, AbstractMatrix},
    PS,
    P
} <: AbstractPoincareInvariant
    Ω::ΩT  # symplectic matrix or function returning one
    pointspec::PS  # specifies how many points
    plan::P  # plan for chebyshev transform, differentiation, etc...
end

function SecondPoincareInvariant{T, D}(
    Ω::ΩT, N, P::Type=DEFAULT_PLAN_TYPE
) where {T, D, ΩT}
    ps = getpointspec(N, P)
    plan = P{T, D}(Ω, ps)
    SecondPoincareInvariant{T, D, ΩT, typeof(ps), typeof(plan)}(Ω, ps, plan)
end

function SecondPoincareInvariant{T, D}(Ω::ΩT, ps::PS, plan::P) where {T, D, ΩT, PS, P}
    SecondPoincareInvariant{T, D, ΩT, PS, P}(Ω, ps, plan)
end

const PI2 = SecondPoincareInvariant

getdim(::SecondPoincareInvariant{<:Any, D, <:Any, <:Any, <:Any}) where D = D
getform(pinv::SecondPoincareInvariant) = pinv.Ω

# Interface should define getpoints(f, T, N, P) method
function getpoints(
    f::Function,
    pinv::SecondPoincareInvariant{T, D, ΩT, PS, P}
) where {T, D, ΩT, PS, P}
    getpoints(f, T, pinv.pointspec, P)
end

getpoints(pinv::SecondPoincareInvariant) = getpoints((x, y) -> (x, y), pinv)

getpointspec(pinv::SecondPoincareInvariant) = pinv.pointspec

## Implementation(s)

include("Chebyshev.jl")
include("FiniteDifferences.jl")

using .Chebyshev: ChebyshevPlan
using .FiniteDifferences: FiniteDiffPlan

const DEFAULT_PLAN_TYPE = ChebyshevPlan

end  # module SecondPoincareInvariants
