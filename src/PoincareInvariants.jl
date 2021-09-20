module PoincareInvariants

using Reexport

export AbstractPoincareInvariant, compute

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

@reexport module FirstPoincareInvariants
    include("FirstPoincareInvariants.jl")
end

@reexport module SecondPoincareInvariants
    include("SecondPoincareInvariants/SecondPoincareInvariants.jl")
end

end  # module
