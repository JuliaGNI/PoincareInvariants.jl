@reexport module FirstPoincareInvariants

import ...PoincareInvariants: getpointspec, getform

using ...PoincareInvariants: AbstractPoincareInvariant

export FirstPoincareInvariant, PI1

struct FirstPoincareInvariant{T, D, θT, P} <: AbstractPoincareInvariant
    θ::θT
    N::Int
    plan::P
end

function FirstPoincareInvariant{T, D}(
    θ::θT, N::Integer, P::Type=DEFAULT_PLAN_TYPE
) where {T, D, θT}
    ps = getpointspec(N, P)
    plan = P{T, D}(θ, N)
    FirstPoincareInvariant{T, D, θT, typeof(plan)}(Ω, ps, plan)
end

FirstPoincareInvariant{T, D}(θ::θT, N::Integer, plan::P) where {T, D, θT, P} =
    FirstPoincareInvariant{T, D, θT, P}(θ, N, plan)

const PI1 = FirstPoincareInvariant

# Interface

getpointspec(pinv::FirstPoincareInvariant) = pinv.N
getform(pinv::FirstPoincareInvariant) = pinv.θ

## Implementations

include("LoopFinDiff.jl")

using .LoopFinDiff: LoopFinDiffPlan

# TODO: rename this to PI1_DEFAULT_PLAN_TYPE
const DEFAULT_PLAN_TYPE = LoopFinDiffPlan

end  # module FirstPoincareInvariants
