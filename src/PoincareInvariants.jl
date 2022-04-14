"""
    PoincareInvariants

A Julia library for the computation of Poincaré integral invariants.
"""
module PoincareInvariants

export AbstractPoincareInvariant, compute!
export FirstPoincareInvariant, CanonicalFirstPI, FirstPI
export SecondPoincareInvariant, CanonicalSecondPI, SecondPI

export getpoints, getpointspec, getpointnum
export getdim, getform, getplan

export FirstFinDiffPlan, FirstFourierPlan
export SecondChebyshevPlan, SecondFinDiffPlan

export canonical_one_form, CanonicalSymplecticMatrix, canonical_two_form

"""
    AbstractPoincareInvariant{T, D}

represents a Poincare integral invariant in a phase space of dimension D using numeric type
T for calculations.
"""
abstract type AbstractPoincareInvariant{T, D} end

"""
    compute!(pinv::AbstractPoincareInvariant, args...)

computes a Poincaré invariant.
"""
function compute! end

"""
    getplan(pinv::AbstractPoincareInvariant)

returns plan that will be used to `compute!` `pinv`.
"""
function getplan end

"""
    getpointspec(pinv::AbstractPoincareInvariant)

get point specification, which may, for example, be a tuple specifying a grid or
a number giving the number of points used to sample in phase space.
"""
function getpointspec end

"""
    getpoints(pinv::AbstractPoincareInvariant)

returns points on which to evaluate the phase space line or surface parameterisation
so as to `compute!` `pinv`.
"""
function getpoints(f::Function, pinv::AbstractPoincareInvariant{T, D}) where {T, D}
    getpoints(f, T, getpointspec(pinv), typeof(getplan(pinv)))
end

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

using .CanonicalSymplecticForms: canonical_one_form, CanonicalSymplecticMatrix,
    canonical_two_form


## FirstPoincareInvariant ##

struct FirstPoincareInvariant{T, D, θT, P} <: AbstractPoincareInvariant{T, D}
    θ::θT
    N::Int
    plan::P
end

function FirstPoincareInvariant{T, D}(
    θ::θT, N::Integer, P::Type=DEFAULT_FIRST_PLAN
) where {T, D, θT}
    ps = getpointspec(N, P)
    plan = P{T, D}(θ, N)
    FirstPoincareInvariant{T, D, θT, typeof(plan)}(θ, ps, plan)
end

FirstPoincareInvariant{T, D}(θ::θT, N::Integer, plan::P) where {T, D, θT, P} =
    FirstPoincareInvariant{T, D, θT, P}(θ, N, plan)

function FirstPoincareInvariant{T, D, typeof(canonical_one_form)}(
    N, P=DEFAULT_FIRST_PLAN
) where {T, D}
    FirstPoincareInvariant{T, D}(canonical_one_form, N, P)
end

const FirstPI = FirstPoincareInvariant
const CanonicalFirstPI{T, D} = FirstPoincareInvariant{T, D, typeof(canonical_one_form)}

# Interface

getpointspec(pinv::FirstPoincareInvariant) = pinv.N
getpoints(pinv::FirstPoincareInvariant) = getpoints(identity, pinv)
getform(pinv::FirstPoincareInvariant) = pinv.θ
getplan(pinv::FirstPoincareInvariant) = pinv.plan
getdim(::FirstPoincareInvariant{<:Any, D}) where D = D

# Implementations

include("FirstFinDiffPlans.jl")
include("FirstFourierPlans.jl")

using .FirstFinDiffPlans: FirstFinDiffPlan
using .FirstFourierPlans: FirstFourierPlan

const DEFAULT_FIRST_PLAN = FirstFourierPlan

## SecondPoincareInvariant ##

struct SecondPoincareInvariant{
    T,  # phase space and return type
    D,  # dimension of phase space
    ωT,
    PS,
    P
} <: AbstractPoincareInvariant{T, D}
    ω::ωT  # symplectic matrix or function returning one
    pointspec::PS  # specifies how many points
    plan::P  # plan for chebyshev transform, differentiation, etc...
end

function SecondPoincareInvariant{T, D}(
    ω::ωT, N, P::Type=DEFAULT_SECOND_PLAN
) where {T, D, ωT}
    ps = getpointspec(N, P)
    plan = P{T, D}(ω, ps)
    SecondPoincareInvariant{T, D, ωT, typeof(ps), typeof(plan)}(ω, ps, plan)
end

function SecondPoincareInvariant{T, D}(ω::ωT, ps::PS, plan::P) where {T, D, ωT, PS, P}
    SecondPoincareInvariant{T, D, ωT, PS, P}(ω, ps, plan)
end

function SecondPoincareInvariant{T, D, CanonicalSymplecticMatrix{T}}(
    N, P=DEFAULT_SECOND_PLAN
) where {T, D}
    ω = CanonicalSymplecticMatrix{T}(D)
    SecondPoincareInvariant{T, D}(ω, N, P)
end

const SecondPI = SecondPoincareInvariant
const CanonicalSecondPI{T, D} = SecondPoincareInvariant{T, D, CanonicalSymplecticMatrix{T}}

# Interface

getdim(::SecondPoincareInvariant{<:Any, D, <:Any, <:Any, <:Any}) where D = D
getform(pinv::SecondPoincareInvariant) = pinv.ω
getplan(pinv::SecondPoincareInvariant) = pinv.plan

getpoints(pinv::SecondPoincareInvariant) = getpoints((x, y) -> (x, y), pinv)

getpointspec(pinv::SecondPoincareInvariant) = pinv.pointspec

## Implementations

include("SecondChebyshevPlans.jl")
include("SecondFinDiffPlans.jl")

using .SecondChebyshevPlans: SecondChebyshevPlan
using .SecondFinDiffPlans: SecondFinDiffPlan

const DEFAULT_SECOND_PLAN = SecondChebyshevPlan

end  # module
