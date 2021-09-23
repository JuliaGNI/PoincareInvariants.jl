module PoincareInvariants

using Reexport

export AbstractPoincareInvariant, compute, getpoints

"""
    AbstractPoincareInvariant

supertype of `FirstPoincareInvariant` and `SecondPoincareInvariant`.
"""
abstract type AbstractPoincareInvariant end

"""
    compute(pinv::AbstractPoincareInvariant, args...)

computes a Poincaré invariant. 
"""
function compute end

"""
    getpoints(pinv::AbstractPoincareInvariant)

returns points on which to evaluate the phase space line or surface parameterisation
so as to `compute` `pinv`.
"""
function getpoints end

"""
    getpointnum(pinv::AbstractPoincareInvariant)

returns number of points to sample in phase space to `compute` `pinv`.
"""
function getpointnum end



@reexport module FirstPoincareInvariants
    include("FirstPoincareInvariants.jl")
end

@reexport module SecondPoincareInvariants
    include("SecondPoincareInvariants/SecondPoincareInvariants.jl")
end

end  # module
