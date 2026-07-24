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

# JuliaGNI ecosystem integration

using GeometricEquations: EnsembleProblem, ODE, HODE, PODE, IODE, LODE,
    equation, timespan, timestep, parameters
using GeometricSolutions: EnsembleSolution, GeometricSolution, nsamples, ntime

export PIEnsembleProblem

include("utils.jl")

## AbstractPoincareInvariant ##

"""
    AbstractPoincareInvariant{T, D}

represents a Poincare integral invariant in a phase space of dimension `D`` using numeric type
`T` for calculations.
"""
abstract type AbstractPoincareInvariant{T, D} end

"""
    compute!(pinv::AbstractPoincareInvariant, points::AbstractMatrix, t::Real=NaN, p=nothing)

computes a Poincaré invariant using setup object `pinv`, where `points` represents the image
of the curve or surface parameterisation evaluated on some set of points, e.g. a grid.

`t` is the time at which the invariant is evaluated `p` are any user supplied optional
arguments. Both `t` and `p` are passed directly to the differential form.

Plan implementations should define a method `compute!(pinv, t::Real, p)`, which acts on the
internal points storage `pinv.points`.
"""
function compute!(
    pinv::AbstractPoincareInvariant,
    points::AbstractMatrix,
    t::Real=NaN, p=nothing
)
    pinv.points .= points
    compute!(pinv, t, p)
end
# use NaN instead of 0.0 for default time so error is thrown in case user provided form
# requires a time. NaN also works in case user required times are Real or FLoat64

_gettime(t::Real, ::Integer)::Real = t
_gettime(ts::AbstractVector{<:Real}, i::Integer)::Real = ts[i]

"""
    compute!(pinv::AbstractPoincareInvariant, points::AbstractVector{<:AbstractVector},
        times::Union{AbstractVector{<:Real}, Real}=NaN, p=nothing)

computes a Poincaré invariant, using the setup object `pinv`, at each time for a set of
trajectories given by `points`.

`points` represents an AbstractVector of trajectories. Each element of `points` is a
trajectory, which is itself an `AbstractVector` of some kind of iterable or vector.
`times` is either a constant, like `NaN`, or a vector of the times at which the
trajectories have been evaluated, i.e. the i-th phase space position in each trajectory
was evaluated at time `times[i]`. `p` is an arbitrary optional parameter which is passed
to the differential form, like the the time.
"""
function compute!(
    pinv::AbstractPoincareInvariant,
    points::AbstractVector{<:AbstractVector},
    times::Union{AbstractVector{<:Real}, Real}=NaN,
    p=nothing
)
    indices = eachindex(points[1])
    @argcheck all(idx -> idx == indices, eachindex.(points)) "indices of vectors of points must match"
    if times isa AbstractVector
        @argcheck indices == eachindex(times) "indices of times and vectors of points must match"
    end
    map(indices) do i
        for (j, v) in enumerate(points)
            pinv.points[j, :] .= v[i]
        end
        compute!(pinv, _gettime(times, i), p)
    end
end

"""
    getplan(pinv::AbstractPoincareInvariant)

returns the plan that will be used to `compute!` `pinv`.
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

"""
    FirstPoincareInvariant{T, D, θT, P} <: AbstractPoincareInvariant{T, D}
    FirstPI{T, D, θT, P} <: AbstractPoincareInvariant{T, D}

setup object to compute the first invariant.
"""
struct FirstPoincareInvariant{T, D, θT, P} <: AbstractPoincareInvariant{T, D}
    "differential one form θ"
    θ::θT
    "N, the number of points to be used"
    N::Int
    "plan used to calculate the invariant"
    plan::P
    "internal point storage"
    points::Matrix{T}
end

"""
    FirstPoincareInvariant{T, D}(θ::θT, N::Integer, plan::P)

constructs a `FirstPoincareInvariant` setup object to calculate integral invariants, given
a numeric type `T`, a phase space dimension `D`, a differential form `θ`, a number of points
`N` and a `plan`. Note that the number of points may not be exactly `N`. The true number
used depends on the implementation and is guaranteed to be no smaller than `N` and not too
much larger.
"""
function FirstPoincareInvariant{T, D}(θ::θT, N::Integer, plan::P) where {T, D, θT, P}
    n = getpointspec(N, P) |> getpointnum
    points = Matrix{T}(undef, n, D)
    FirstPoincareInvariant{T, D, θT, P}(θ, n, plan, points)
end

"""
    FirstPoincareInvariant{T, D}(θ::θT, N::Integer, P::Type=DEFAULT_FIRST_PLAN)

constructs a `FirstPoincareInvariant` setup object to calculate integral invariants, given
a numeric type `T`, a phase space dimension `D`, a differential form `θ`, a number of points
`N` and a `plan` type `P`, to be initialised and used for future computation.

The plan type defaults to `DEFAULT_FIRST_PLAN`, which is currently set to `FirstFourierPlan`.
"""
function FirstPoincareInvariant{T, D}(
    θ::θT, N::Integer, P::Type=DEFAULT_FIRST_PLAN
) where {T, D, θT}
    plan = P{T, D}(θ, getpointspec(N, P))
    FirstPoincareInvariant{T, D}(θ, N, plan)
end

"""
    FirstPoincareInvariant{T, D, typeof(canonical_one_form)}(N::Integer, P=DEFAULT_FIRST_PLAN)
    CanonicalFirstPI{T, D}(N::Integer, P=DEFAULT_FIRST_PLAN)

creates a setup object to compute the first integral invariant using the canonical one form
in phase space of dimension `D` with numeric type `T`, `N` points and plan type `P`.
"""
function FirstPoincareInvariant{T, D, typeof(canonical_one_form)}(
    N::Integer, P=DEFAULT_FIRST_PLAN
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

"""
    SecondPoincareInvariant{T, D, ωT, PS, P} <: AbstractPoincareInvariant{T, D}
    SecondPI{T, D, ωT, PS, P} <: AbstractPoincareInvariant{T, D}

setup object used to compute the second invariant.
"""
struct SecondPoincareInvariant{
    T,  # phase space and return type
    D,  # dimension of phase space
    ωT,
    PS,
    P
} <: AbstractPoincareInvariant{T, D}
    "differential two form to be integrated"
    ω::ωT
    "specification of points at which to sample surface parameterisation"
    pointspec::PS
    "plan to use for computation of invariant"
    plan::P
    "internal point storage"
    points::Matrix{T}
end

"""
    SecondPoincareInvariant{T, D}(ω, N, plan::P)

constructs a `SecondPoincareInvariant` setup object to calculate integral invariants, given
a numeric type `T`, a phase space dimension `D`, a differential form `ω`, a point
specification `N`, usually a number of points, and a `plan`.
"""
function SecondPoincareInvariant{T, D}(ω::ωT, N, plan::P) where {T, D, ωT, P}
    ps = getpointspec(N, P)
    points = Matrix{T}(undef, getpointnum(ps), D)
    SecondPoincareInvariant{T, D, ωT, typeof(ps), P}(ω, ps, plan, points)
end

"""
    SecondPoincareInvariant{T, D}(ω, N, P::Type=DEFAULT_SECOND_PLAN)

constructs a `SecondPoincareInvariant` setup object to calculate integral invariants, given
a numeric type `T`, a phase space dimension `D`, a differential form `ω`, a point
specification `N` and a `plan` type `P`, to be initialised and used for future computation.

The plan type defaults to `DEFAULT_SECOND_PLAN`, which is currently set to
`SecondChebyshevPlan`. Note that the number of points may not be exactly `N`. The true
number used depends on the implementation and is guaranteed to be no smaller than `N` and
not too much larger. For the `SecondFinDiffPlan`, the grid of points may also be specified
as a tuple `(Nx, Ny)`, if non-square grids are sought.
"""
function SecondPoincareInvariant{T, D}(
    ω, N, P::Type=DEFAULT_SECOND_PLAN
) where {T, D}
    plan = P{T, D}(ω, getpointspec(N, P))
    SecondPoincareInvariant{T, D}(ω, N, plan)
end

"""
    SecondPoincareInvariant{T, D, CanonicalSymplecticMatrix{T}}(N, P=DEFAULT_FIRST_PLAN)
    CanonicalSecondPI{T, D}(N, P=DEFAULT_FIRST_PLAN)

creates a setup object to compute the second integral invariant using the canonical two form
in phase space of dimension `D` with numeric type `T`, point specification `N` and plan
type `P`.
"""
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

## Integration ##

"""
    PIEnsembleProblem(prob, pinv, init)

builds a `GeometricEquations.EnsembleProblem` from a base problem `prob`, a Poincaré
invariant setup object `pinv`, and a curve or surface parameterisation `init`.

`prob` is a geometric `EquationProblem` (e.g. an `ODEProblem`, `HODEProblem`, `IODEProblem`
or `LODEProblem`) specifying the dynamics, time span and time step. The parameterisation
`init` is sampled with [`getpoints`](@ref), and each sampled phase space point becomes the
initial condition of one member of the ensemble. For `HODEProblem`/`PODEProblem` the point
is split into the position and momentum halves `(q, p)`; for `IODEProblem`/`LODEProblem` the
momentum is initialised from the equation's one-form as `p₀ = ϑ(t₀, q₀)`.

The returned ensemble can be integrated with `GeometricIntegrators.integrate` and the
invariant computed from the resulting solution with [`compute!`](@ref).
"""
function PIEnsembleProblem(prob, pinv::AbstractPoincareInvariant, init)
    points = getpoints(init, pinv)
    equ = equation(prob)
    t₀ = timespan(prob)[begin]
    par = parameters(prob)
    ics = [_ensemble_ics(equ, view(points, i, :), t₀, par) for i in axes(points, 1)]
    EnsembleProblem(equ, timespan(prob), timestep(prob), ics, par)
end

# map a single D-dimensional phase space point to the initial condition NamedTuple
# expected by the underlying geometric equation
_ensemble_ics(::ODE, point, t₀, params) = (q = collect(point),)

function _ensemble_ics(::Union{HODE, PODE}, point, t₀, params)
    mid = length(point) ÷ 2
    (q = collect(point[1:mid]), p = collect(point[mid+1:end]))
end

function _ensemble_ics(equ::Union{IODE, LODE}, point, t₀, params)
    q = collect(point)
    p = zero(q)
    v = zero(q)
    equ.ϑ(p, t₀, q, v, params)
    (q = q, p = p)
end

"""
    compute!(pinv::AbstractPoincareInvariant, sol::EnsembleSolution, p=nothing)

computes the Poincaré invariant at every saved time step of an `EnsembleSolution` returned
by `GeometricIntegrators.integrate`. Each ensemble member provides the trajectory of one
point on the curve or surface. Returns a `Vector` holding one invariant value per time step.

The optional parameter `p` is passed to the differential form, just like the time.
"""
function compute!(
    pinv::AbstractPoincareInvariant,
    sol::EnsembleSolution,
    p=nothing
)
    D = getdim(pinv)
    n = nsamples(sol)
    # assumes the time grid is the same for all ensemble members
    m1 = sol[1]
    hasp = hasproperty(m1, :p)
    map(0:ntime(m1)) do i
        for j in 1:n
            m = sol[j]
            q = m.q[i]
            if hasp && length(q) < D
                # canonical (q, p) phase space: stack position and momentum halves
                mid = length(q)
                @views pinv.points[j, 1:mid] .= q
                @views pinv.points[j, mid+1:D] .= m.p[i]
            else
                # phase space coincides with the position variable q
                @views pinv.points[j, 1:D] .= q
            end
        end
        compute!(pinv, m1.t[i], p)
    end
end

end  # module
