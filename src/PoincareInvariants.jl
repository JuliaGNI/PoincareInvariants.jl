module PoincareInvariants

# For copy of paduatransform!
using LinearAlgebra: rmul!
using FastTransforms: paduavalsmat, trianglecfsvec!

# For the Padua transforms and manipulations in Chebyshev base
using ApproxFunOrthogonalPolynomials
using FastTransforms: PaduaTransformPlan, plan_paduatransform!, paduapoints

using LinearAlgebra: dot

# to represent points
using StaticArrays

# Callable is Union{Function, Type}
using Base: Callable

export get_padua_points
export AbstractPoincareInvariant, PoincareInvariant1, PoincareInvariant2, compute

include("padua.jl")

abstract type AbstractPoincareInvariant end

include("PoincareInvariant2.jl")

end  # module
