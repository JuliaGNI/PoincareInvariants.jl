@reexport module SecondPoincareInvariants

# Callable is Union{Function, Type}
using Base: Callable

using ..PoincareInvariants: AbstractPoincareInvariant, @argcheck
import ..PoincareInvariants: compute!, getpoints, getpointnum, getdim

export SecondPoincareInvariant

## Implementations

include("ChebyshevImplementation.jl")

## SecondPoincareInvariant ##

# TODO: add isinplace check

struct SecondPoincareInvariant{
	T,  # phase space and return type
	ΩT <: Union{Callable, AbstractMatrix}, 
	S
} <: AbstractPoincareInvariant
	Ω::ΩT  # symplectic matrix or function returning one
	D::Int
	setup::S  # setup object for preallocating operators, transform plans, etc.
end

function SecondPoincareInvariant{T}(Ω::ΩT, D::Integer, N::Integer) where {T, ΩT <: AbstractMatrix}
	@argcheck size(Ω) == (D, D) "Ω must be a $D × $D matrix"
	setup = ChebyshevImplementation.ChebyshevSetup{T}(Ω, D, N)
	SecondPoincareInvariant{T, ΩT, typeof(setup)}(Ω, D, setup)
end

function SecondPoincareInvariant{T}(
	Ω::ΩT, D::Integer, N::Integer, ::Val{inplace}
) where {T, ΩT <: Callable, inplace}
	setup = ChebyshevImplementation.ChebyshevSetup{T}(Ω, D, N, Val(inplace))
	SecondPoincareInvariant{T, ΩT, typeof(setup)}(Ω, D, setup)
end

# Unexported convenience alias
const PI2 = SecondPoincareInvariant

getdim(pinv::SecondPoincareInvariant) = pinv.D

## Internal interface to implementation ##
compute!(pinv::SecondPoincareInvariant, args...) = _compute(pinv.setup, pinv.Ω, pinv.D, args...)
getpoints(pinv::SecondPoincareInvariant) = getpoints(pinv.setup)
getpointnum(pinv::SecondPoincareInvariant) = getpointnum(pinv.setup)

end  # module SecondPoincareInvariants
