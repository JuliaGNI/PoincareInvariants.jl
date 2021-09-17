module PoincareInvariants

using Reexport

export AbstractPoincareInvariant, compute

"""
    AbstractPoincareInvariant

supertype of setup objects `FirstPoincareInvariant` and `SecondPoincareInvariant`.
"""
abstract type AbstractPoincareInvariant end

"""
    compute(pinv::AbstractPoincareInvariant, args...)

computes a Poincaré invariant.
"""
function compute end

@reexport module CanonicalSymplecticMatrices include("CanonicalSymplecticMatrices.jl") end
@reexport module FirstPoincareInvariants include("FirstPoincareInvariants.jl") end
@reexport module SecondPoincareInvariants include("SecondPoincareInvariants.jl") end

end  # module
