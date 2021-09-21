using ApproxFunOrthogonalPolynomials
using FastTransforms: PaduaTransformPlan, plan_paduatransform!
using ApproxFunBase: TransformPlan, ITransformPlan, plan_transform, plan_itransform

using BlockBandedMatrices: BandedBlockBandedMatrix
using BlockArrays: BlockRange
using StaticArrays: SVector
using SparseArrays: sparse

using LinearAlgebra: mul!, dot

# Callable is Union{Function, Type}
using Base: Callable

using ..PoincareInvariants: AbstractPoincareInvariant
import ..PoincareInvariants: compute, getpoints

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
	getsetup(D::Integer, T::Type{T}, Ω::Union{Callable, AbstractMatrix}, N::Integer)

gets setup object for preallocating arrays, transform plans and operators.
"""
function getsetup end

struct SecondPoincareInvariant{
	D, T,  # phase space dimension and type
	ΩT <: Union{Callable, AbstractMatrix}, 
	S
} <: AbstractPoincareInvariant
	Ω::ΩT  # symplectic matrix
	N::Int  # number of points
	setup::S  # setup object for preallocating operators, transform plans, etc.

	function SecondPoincareInvariant{D, T}(
		Ω::ΩT, N::Integer, setup::S
	) where {D, T, ΩT <: Union{Callable, AbstractMatrix}, S}
		checkΩ(Ω, D, T)
		checkpaduanum(N)
		new{D, T, ΩT, S}(Ω, N, setup)
	end
end

function SecondPoincareInvariant{D, T}(
	Ω::ΩT, N::Integer
) where {D, T, ΩT <: Union{Callable, AbstractMatrix}}
	paduanum = nextpaduanum(N)
	SecondPoincareInvariant{D, T}(Ω, paduanum, getsetup(D, T, Ω, paduanum))
end

const PI2 = SecondPoincareInvariant

include("Chebyshev.jl")

getpoints(pinv::SecondPoincareInvariant{<:Any, T}) where T = getpaduapoints(T, ceil(Int, getdegree(pinv.N)))

include("implementation1.jl")
include("implementation2.jl")

getsetup(D::Integer, ::Type{T}, Ω::AbstractMatrix, N::Integer) where T = Setup2(D, T, Ω, N)

include("CanonicalSymplecticMatrices.jl")
