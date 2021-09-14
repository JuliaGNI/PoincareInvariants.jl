module PoincareInvariants

# For copy of paduatransform!
using LinearAlgebra: rmul!
using FastTransforms: paduavalsmat, trianglecfsvec!

# For the Padua transforms and manipulations in Chebyshev base
using ApproxFunBase: TransformPlan, ITransformPlan, plan_transform, plan_itransform
using ApproxFunOrthogonalPolynomials
using FastTransforms: PaduaTransformPlan, plan_paduatransform!, paduapoints

using BlockArrays: BlockRange

# For operator preallocation
using ApproxFunBase: BandedBlockBandedMatrix, FiniteRange

using LinearAlgebra

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
