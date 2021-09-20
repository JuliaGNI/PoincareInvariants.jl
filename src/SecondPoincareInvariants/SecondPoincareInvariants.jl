using ApproxFunOrthogonalPolynomials
using FastTransforms: PaduaTransformPlan, plan_paduatransform!
using ApproxFunBase: TransformPlan, ITransformPlan, plan_transform, plan_itransform

using BlockBandedMatrices: BandedBlockBandedMatrix
using BlockArrays: BlockRange
using StaticArrays: SVector

using LinearAlgebra: mul!, dot

# Callable is Union{Function, Type}
using Base: Callable

using ..PoincareInvariants: AbstractPoincareInvariant
import ..PoincareInvariants: compute

export SecondPoincareInvariant
export getpoints, nextpointnum

function checkΩ(Ω::AbstractMatrix, D::Integer, ::Type)
	size(Ω) == (D, D) || throw(ArgumentError(
		"Ω must be a $D × $D AbstractMatrix, not $(size(Ω, 1)) × $(size(Ω, 2))"
	))
end

function checkΩ(Ω::Callable, D::Integer, ::Type{T}) where T
	ω = Ω(zeros(T, D))
	ω isa AbstractMatrix || throw(ArgumentError("Ω must return some AbstractMatrix"))
	checkΩ(ω, D, T)
end

"""
	makesetup(D::Integer, T::Type{T}, N::Integer, Ω::Union{Callable, AbstractMatrix})

creates setup object for preallocating arrays, transform plans and operators.
"""
function makesetup end

struct SecondPoincareInvariant{
	D, T,  # phase space dimension and type
	ΩT <: Union{Callable, AbstractMatrix}, 
	S} <: AbstractPoincareInvariant
	N::Int  # number of points
	Ω::ΩT  # symplectic matrix
	setup::S  # setup object for preallocating operators, transform plans, etc.
end

# The user should use this method
function SecondPoincareInvariant{D, T}(N::Integer, Ω::ΩT) where {D, T, ΩT}
	checkΩ(Ω)

	setup = makesetup(D, T, N, Ω)
	S = typeof(setup)

	SecondPoincareInvariant{D, T, ΩT, S}(N, Ω, setup)
end

const PI2 = SecondPoincareInvariant

include("Chebyshev.jl")

# include("implementation1.jl")

makesetup(D::Integer, ::Type{T}, N::Integer, Ω::AbstractMatrix) where T = Setup1(D, T, N, Ω)

include("CanonicalSymplecticMatrices.jl")
