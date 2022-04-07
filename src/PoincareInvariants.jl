"""
    PoincareInvariants

A Julia library for the computation of Poincaré integral invariants.
"""
module PoincareInvariants

# Callable is Union{Function, Type}
using Base: Callable

export AbstractPoincareInvariant, compute!
export FirstPoincareInvariant, PI1, FirstFinDiffPlan
export SecondPoincareInvariant, PI2, SecondChebyshevPlan, SecondFinDiffPlan

export getpoints, getpointspec, getpointnum
export getdim, getform

export CanonicalSymplecticTwoForm

"""
    AbstractPoincareInvariant

supertype of `FirstPoincareInvariant` and `SecondPoincareInvariant`.
"""
abstract type AbstractPoincareInvariant end

"""
    compute!(pinv::AbstractPoincareInvariant, args...)

computes a Poincaré invariant.
"""
function compute! end

"""
    getpoints(pinv::AbstractPoincareInvariant)

returns points on which to evaluate the phase space line or surface parameterisation
so as to `compute!` `pinv`.
"""
function getpoints end

"""
    getpointspec(pinv::AbstractPoincareInvariant)

get point specification, which may, for example, be a tuple specifying a grid or
a number giving the number of points used to sample in phase space.
"""
function getpointspec end

"""
    getpointnum(pinv::AbstractPoincareInvariant)

returns number of points to sample in phase space to `compute!` `pinv`.
"""
getpointnum(pinv::AbstractPoincareInvariant) = getpointnum(getpointspec(pinv))
getpointnum(N::Integer) = N

"""
    getdim(pinv::AbstractPoincareInvariant)

returns dimension of phase space to `compute!` `pinv` in.
"""
function getdim end

"""
    getform(pinv::AbstractPoincareInvariant)

get invariant one- or two-form.
"""
function getform end


## Utils ##

include("utils.jl")
include("CanonicalSymplecticForms.jl")

using .CanonicalSymplecticForms: CanonicalSymplecticTwoForm


## FirstPoincareInvariant ##

struct FirstPoincareInvariant{T, D, θT, P} <: AbstractPoincareInvariant
    θ::θT
    N::Int
    plan::P
end

function FirstPoincareInvariant{T, D}(
    θ::θT, N::Integer, P::Type=DEFAULT_FIRST_PLAN
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

# Implementations

include("FirstFinDiffPlans.jl")

using .FirstFinDiffPlans: FirstFinDiffPlan

const DEFAULT_FIRST_PLAN = FirstFinDiffPlan

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
    Ω::ΩT, N, P::Type=DEFAULT_SECOND_PLAN
) where {T, D, ΩT}
    ps = getpointspec(N, P)
    plan = P{T, D}(Ω, ps)
    SecondPoincareInvariant{T, D, ΩT, typeof(ps), typeof(plan)}(Ω, ps, plan)
end

function SecondPoincareInvariant{T, D}(Ω::ΩT, ps::PS, plan::P) where {T, D, ΩT, PS, P}
    SecondPoincareInvariant{T, D, ΩT, PS, P}(Ω, ps, plan)
end

const PI2 = SecondPoincareInvariant

# Interface

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

## Implementations

include("SecondChebyshevPlans.jl")
include("SecondFinDiffPlans.jl")

using .SecondChebyshevPlans: SecondChebyshevPlan
using .SecondFinDiffPlans: SecondFinDiffPlan

const DEFAULT_SECOND_PLAN = SecondChebyshevPlan

end  # module
