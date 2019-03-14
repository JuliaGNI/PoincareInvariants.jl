
function vector_matrix_vector_product(x::Vector{T}, A::Matrix{T}, y::Vector{T}) where {T}
@assert length(x) == size(A,1)
@assert length(y) == size(A,2)

local b::T = 0
local c::T = 0

@inbounds for i in 1:length(x)
    c = 0
    for j in 1:length(y)
        c += A[i,j] * y[j]
    end
    b += c * x[i]
end

return b
end

function vector_matrix_vector_product(x::Matrix{T}, A::Matrix{T}, y::Matrix{T}, ind::Int) where {T}
@assert size(x,1) == size(A,1)
@assert size(y,1) == size(A,2)
@assert ind ≤ size(x,2)
@assert ind ≤ size(y,2)

local b::T = 0
local c::T = 0

@inbounds for i in 1:size(x,1)
    c = 0
    for j in 1:size(y,1)
        c += A[i,j] * y[j,ind]
    end
    b += c * x[i,ind]
end

return b
end

function vector_vector_product(x::Matrix{T}, y::Matrix{T}, ind::Int) where {T}
@assert size(x,1) == size(y,1)
@assert ind ≤ size(x,2)
@assert ind ≤ size(y,2)

local b::T = 0

@inbounds for i in 1:size(x,1)
    b += x[i,ind] * y[i,ind]
end

return b
end
