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

include("utils.jl")

## AbstractPoincareInvariant ##

"""
    AbstractPoincareInvariant{T, D}

represents a Poincare integral invariant in a phase space of dimension D using numeric type
T for calculations.
"""
abstract type AbstractPoincareInvariant{T, D} end

# these types will be accepted as a point cloud representing an area or line
# Any other type will get fed to a PoincareInvariantIter constructor
# The rows of an AbstractMatrix are points and
# the elements of the AbstractVector{AbstractVector} are points
const PointCollection{T} = Union{AbstractMatrix{T}, AbstractVector{AbstractVector{T}}}

writepoints!(pinv::AbstractPoincareInvariant, points::AbstractMatrix) = (pinv.points .= points)
function writepoints!(pinv::AbstractPoincareInvariant, points::AbstractVector{AbstractVector})
    for i in 1:length(points)
        pinv.points[:, i] .= points[i]
    end
end

"""
    compute!(pinv::AbstractPoincareInvariant, args...)

computes a Poincaré invariant.

implementations should define a method `compute!(pinv, t, p)`, which acts on the internal
points storage.
"""
function compute!(
    pinv::AbstractPoincareInvariant,
    points::PointCollection,
    t::Real, p
)
    writepoints!(pinv, points)
    compute!(pinv, t, p)
end

function compute!(pinv::AbstractPoincareInvariant, data, times, p)
    @argcheck eltype(data) <: PointCollection
    map(enumerate(data)) do (i, points)
        writepoints!(pinv, points)
        compute!(pinv, times[i], p)
    end
end

#=
function compute!(picache::ThreadCache{<:AbstractPoincareInvariant{T}}, data, ts, p) where T
    n = datalength(data)
    lk = ReentrantLock()

    out = Vector{T}(undef, n)
    @threads for i in 1:n
        # each thread gets one pinv
        pinv = picache[Threads.threadid()]

        # lock in case data has mutable state
        lock(_ -> writepoints!(pinv, data, i), lk)

        out[i] = compute!(pinv, ts[i], p)
    end

    out
end

function compute!(pinv::AbstractPoincareInvariant, data, t, p)
    iter = PoincareInvariantIter(pinv, pointsiter, args...)
    map(iter) do (points, t, p)
        compute!(pinv, points, t, p)
    end
end
=#

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


## Canonical Forms ##

include("CanonicalSymplecticForms.jl")

using .CanonicalSymplecticForms: canonical_one_form, CanonicalSymplecticMatrix,
    canonical_two_form


## FirstPoincareInvariant ##

struct FirstPoincareInvariant{T, D, θT, P} <: AbstractPoincareInvariant{T, D}
    θ::θT
    N::Int
    plan::P
    points::Matrix{T}
end

function FirstPoincareInvariant{T, D}(θ::θT, N::Integer, plan::P) where {T, D, θT, P}
    n = getpointspec(N, P) |> getpointnum
    points = Matrix{T}(undef, n, D)
    FirstPoincareInvariant{T, D, θT, P}(θ, n, plan, points)
end

function FirstPoincareInvariant{T, D}(
    θ::θT, N::Integer, P::Type=DEFAULT_FIRST_PLAN
) where {T, D, θT}
    plan = P{T, D}(θ, getpointspec(N, P))
    FirstPoincareInvariant{T, D}(θ, N, plan)
end

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
    points::Matrix{T}
end

function SecondPoincareInvariant{T, D}(ω::ωT, N, plan::P) where {T, D, ωT, P}
    ps = getpointspec(N, P)
    points = Matrix{T}(undef, getpointnum(ps), D)
    SecondPoincareInvariant{T, D, ωT, typeof(ps), P}(ω, ps, plan, points)
end

function SecondPoincareInvariant{T, D}(
    ω, N, P::Type=DEFAULT_SECOND_PLAN
) where {T, D}
    plan = P{T, D}(ω, getpointspec(N, P))
    SecondPoincareInvariant{T, D}(ω, N, plan)
end

function SecondPoincareInvariant{T, D, CanonicalSymplecticMatrix{T}}(
    N, P::Type=DEFAULT_SECOND_PLAN
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

# Implementations

include("SecondChebyshevPlans.jl")
include("SecondFinDiffPlans.jl")

using .SecondChebyshevPlans: SecondChebyshevPlan
using .SecondFinDiffPlans: SecondFinDiffPlan

const DEFAULT_SECOND_PLAN = SecondChebyshevPlan

end  # module
