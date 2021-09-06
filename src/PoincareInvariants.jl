module PoincareInvariants

# For copy of paduatransform!
using LinearAlgebra: rmul!
using FastTransforms: paduavalsmat, trianglecfsvec!

# For the Padua transforms and manipulations in Chebyshev base
using ApproxFunOrthogonalPolynomials: Chebyshev, Fun
using FastTransforms: PaduaTransformPlan, plan_paduatransform!, paduapoints

# to represent points
using StaticArrays

# Callable is Union{Function, Type}
using Base: Callable

export get_padua_points
export AbstractPoincareInvariant, PoincareInvariant1, PoincareInvariant2, calculate

include("padua.jl")

# Do we want PoincareInvariant2 <: AbstractPoincareInvariant2 <: AbstractPoincareInvariant ?
# I guess AbstractPoincareInvariant2 would be worth it if we had different types of
# PoincareInvariant2, like CanonicalPoincareInvariant2 may be
abstract type AbstractPoincareInvariant end

include("PoincareInvariant2.jl")

end  # module
