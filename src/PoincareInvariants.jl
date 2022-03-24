"""
    PoincareInvariants

A Julia library for the computation of Poincaré integral invariants.
"""
module PoincareInvariants

using Reexport

export AbstractPoincareInvariant, compute!, getdim, getform,
    getpoints, getpointspec, getpointnum

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
getpointspec(ps, P) = getpointnum(ps)
# falls back to giving number of points

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

include("utils.jl")

include("CanonicalSymplecticStructures.jl")
include("FirstPoincareInvariants/FirstPoincareInvariants.jl")
include("SecondPoincareInvariants/SecondPoincareInvariants.jl")

end  # module
