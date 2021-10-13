module PoincareInvariants

using Reexport

export AbstractPoincareInvariant, compute!, getdim, getpoints, getpointnum

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
    getpointnum(pinv::AbstractPoincareInvariant)

returns number of points to sample in phase space to `compute!` `pinv`.
"""
function getpointnum end

"""
    getdim(pinv::AbstractPoincareInvariant)

returns dimension of phase space to `compute!` `pinv` in.
"""
function getdim end

include("utils.jl")

include("CanonicalSymplecticStructures.jl")
include("FirstPoincareInvariants.jl")
include("SecondPoincareInvariants/SecondPoincareInvariants.jl")

end  # module
