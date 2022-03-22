module FiniteDifferences

using ...PoincareInvariants: @argcheck
import ...PoincareInvariants: compute!, getpoints, getpointnum
import ..SecondPoincareInvariants: getpointspec

struct FiniteDiffPlan end

function getpointnum(dims::NTuple{2, Integer}, ::Type{<:FiniteDiffPlan})
    nx, ny = getpointspec(dims, FiniteDiffPlan)
    return nx * ny
end

getpointnum(N, ::Type{T}) where T <: FiniteDiffPlan = getpointnum(getpointspec(N, T), T)

nextodd(n::Integer)::Int = n <= 3 ? 3 : 2 * fld(n, 2) + 1

function getpointspec(dims::NTuple{2, Integer}, ::Type{<:FiniteDiffPlan})
    return (nextodd(dims[1]), nextodd(dims[2]))
end

function getpointspec(N::Integer, ::Type{<:FiniteDiffPlan})
    n = nextodd(ceil(Int, sqrt(Int(N))))
    return (n, n)
end

function getpoints(f, ::Type{T}, dims::NTuple{2, Integer}, ::Type{<:FiniteDiffPlan}) where T
    D = length(f(zero(T), zero(T)))
    nx, ny = dims
    N = nx * ny
    out = ntuple(_ -> Vector{T}(undef, N), D)

    i = 1
    for x in range(0, 1, length=nx)
        for y in range(0, 1, length=ny)
            fpnt = f(x, y)
            for d in 1:D
                out[d][i] = fpnt[d]
            end
            i += 1
        end
    end

    # Should return vector [...] instead of ([...],) in 1D case
    return D == 1 ? out[1] : out
end

## Differentiation ##

# convert a cartesian index (x, y) to a linear index and index into vector
_tolin(x::Integer, y::Integer, ny::Integer) = (x - 1) * ny + y

# evaluate derivative of quadratic approxmiation around (x0, y0) at (x0 + x, y0 + y)
# We have approximation f(x, y) = c + cx*x + cy*y + cxx*x^2 + cxy*x*y + cyy*y^2
function _diff(
    vals::AbstractVector{T}, x0::Integer, y0::Integer,
    x::Integer, y::Integer, nx::Integer, ny::Integer,
    ::Val{edge}
) where {T, edge}
    f(ix, iy) = @inbounds vals[_tolin(x0 + ix, y0 + iy, ny)]

    # All coefficients are multiplied by 2*Δx or 2*Δy
    bx = f(1, 0); ax = f(-1, 0)
    fx = cx = bx - ax

    by = f(0, 1); ay = f(0, -1)
    fy = cy = by - ay

    if edge  # meaning if (x, y) not (0, 0)
        f00 = f(0, 0)
        cxx = bx + ax - 2*f00
        cyy = by + ay - 2*f00
        cxy = (f(1, 1) + f(-1, -1) - f(1, -1) - f(-1, 1)) / 2

        fx += 2*cxx * x + cxy * y
        fy += 2*cyy * y + cxy * x
    end

    invΔx = nx - 1; invΔy = ny - 1
    return fx * invΔx / 2, fy * invΔy / 2
end

function differentiate!(
    ∂x::AbstractVector, ∂y::AbstractVector, vals::AbstractVector{T},
    dims::NTuple{2, Integer}
) where T
    nx, ny = dims
    Base.require_one_based_indexing(∂x, ∂y, vals)
    @argcheck nx >= 3 && ny >= 3 "must have at least 3×3 points to calculate both derivatives"
    @argcheck axes(∂x, 1) == axes(∂y, 1) == axes(vals, 1) "axes of derivatives and " *
        "values must match"
    @argcheck length(vals) == nx * ny "length of vals and length implied by dims must match"

    _li(x, y) = _tolin(x, y, ny)

    # corners
    ∂x[_li(1 ,  1)], ∂y[_li(1 ,  1)] = _diff(vals,    2,    2, -1, -1, nx, ny, Val(true))
    ∂x[_li(nx,  1)], ∂y[_li(nx,  1)] = _diff(vals, nx-1,    2, +1, -1, nx, ny, Val(true))
    ∂x[_li(nx, ny)], ∂y[_li(nx, ny)] = _diff(vals, nx-1, ny-1, +1, +1, nx, ny, Val(true))
    ∂x[_li(1 , ny)], ∂y[_li(1 , ny)] = _diff(vals,    2, ny-1, -1, +1, nx, ny, Val(true))

    # edges
    @inbounds for y in 2:ny-1
        ∂x[_li(1 ,  y)], ∂y[_li(1 ,  y)] = _diff(vals,    2, y, -1, 0, nx, ny, Val(true))
        ∂x[_li(nx,  y)], ∂y[_li(nx,  y)] = _diff(vals, nx-1, y, +1, 0, nx, ny, Val(true))
    end

    @inbounds for x in 2:nx-1
        ∂x[_li(x,  1)], ∂y[_li(x,  1)] = _diff(vals, x,    2, 0, -1, nx, ny, Val(true))
        ∂x[_li(x, ny)], ∂y[_li(x, ny)] = _diff(vals, x, ny-1, 0, +1, nx, ny, Val(true))
    end

    # interior (uses central differences)
    @inbounds for x in 2:nx-1
        for y in 2:ny-1
            ∂x[_li(x, y)], ∂y[_li(x, y)] = _diff(vals, x, y, 0, 0, nx, ny, Val(false))
        end
    end

    return ∂x, ∂y
end

function differentiate!(∂x, ∂y, vals, dims::NTuple{2, Integer})
    @argcheck length(∂x) == length(∂y) == length(vals) "number of derivative and vals " *
        "vectors must match"

    for (∂xi, ∂yi, valsi) in zip(∂x, ∂y, vals)
        differentiate!(∂xi, ∂yi, valsi, dims)
    end

    return ∂x, ∂y
end



#=

struct PoincareInvariant2ndTrapezoidal{DT,TT,ET,ΩT} <: AbstractPoincareInvariant2nd{DT}
    equ::ET
    ω::ΩT
    Δt::TT
    nx::Int
    ny::Int
    ntime::Int
    nsave::Int
    nt::Int
    I::OffsetArray{Double64,1,Vector{Double64}}
    J::OffsetArray{Double64,1,Vector{Double64}}
    ΔI::OffsetArray{Double64,1,Vector{Double64}}
    ΔJ::OffsetArray{Double64,1,Vector{Double64}}
end

function PoincareInvariant2ndTrapezoidal(f_equ::Function, f_surface::Function, ω::ΩT, Δt::TT, d::Int, nx::Int, ny::Int, ntime::Int, nsave::Int=1, DT=Float64) where {TT,ΩT}

    if get_config(:verbosity) > 1
        println()
        println("Second Euler-Poincaré Integral Invariant (Trapezoidal)")
        println("======================================================")
        println()
        println(" nx    = ", nx)
        println(" ny    = ", ny)
        println(" ntime = ", ntime)
        println(" nsave = ", nsave)
        println(" Δt    = ", Δt)
        println()
    end

    # compute initial conditions
    q₀ = reshape([f_surface(i/nx, j/ny) for i in 1:nx, j in 1:ny], nx*ny)

    equ = f_equ(q₀)

    # create arrays for results
    nt = div(ntime, nsave)

    I  = OffsetArray(zeros(Double64, nt+1), 0:nt)
    J  = OffsetArray(zeros(Double64, nt+1), 0:nt)
    ΔI = OffsetArray(zeros(Double64, nt+1), 0:nt)
    ΔJ = OffsetArray(zeros(Double64, nt+1), 0:nt)

    PoincareInvariant2ndTrapezoidal{DT,TT,typeof(equ),ΩT}(equ, ω, Δt, nx, ny, ntime, nsave, nt, I, J, ΔI, ΔJ)
end


function interpolate_trajectory(x, i1, j1, λ, μ, γ, nx, ny)
    @assert length(γ) == size(x, 1)
    @assert λ ≥ 0.
    @assert λ ≤ 1.
    @assert μ ≥ 0.
    @assert μ ≤ 1.

    @assert i1 > 0
    @assert i1 < nx
    @assert j1 > 0
    @assert j1 < ny

    i2 = i1 + 1
    j2 = j1 + 1

    for k in eachindex(γ)
        γ[k] = x[k, nx*(j1-1)+i1] * (1-λ) * (1-μ) +
               x[k, nx*(j1-1)+i2] *    λ  * (1-μ) +
               x[k, nx*(j2-1)+i1] * (1-λ) *    μ  +
               x[k, nx*(j2-1)+i2] *    λ  *    μ
    end
end


function interpolate_derivative_i(x, i1, j1, λ, μ, γ̇, nx, ny)
    @assert length(γ̇) == size(x, 1)
    @assert λ ≥ 0.
    @assert λ ≤ 1.
    @assert μ ≥ 0.
    @assert μ ≤ 1.

    @assert i1 > 0
    @assert i1 < nx
    @assert j1 > 0
    @assert j1 < ny

    i2 = i1 + 1
    j2 = j1 + 1

    for k in eachindex(γ̇)
        γ̇[k] = (x[k, nx*(j1-1)+i2] - x[k, nx*(j1-1)+i1]) * (1-μ) +
               (x[k, nx*(j2-1)+i2] - x[k, nx*(j2-1)+i1]) *    μ
    end
end


function interpolate_derivative_j(x, i1, j1, λ, μ, γ̇, nx, ny)
    @assert length(γ̇) == size(x, 1)
    @assert λ ≥ 0.
    @assert λ ≤ 1.
    @assert μ ≥ 0.
    @assert μ ≤ 1.

    @assert i1 > 0
    @assert i1 < nx
    @assert j1 > 0
    @assert j1 < ny

    i2 = i1 + 1
    j2 = j1 + 1

    for k in eachindex(γ̇)
        γ̇[k] = (x[k, nx*(j2-1)+i1] - x[k, nx*(j1-1)+i1]) * (1-λ) +
               (x[k, nx*(j2-1)+i2] - x[k, nx*(j1-1)+i2]) *    λ
    end
end


function integrate(t, γ, γ̇ᵢ, γ̇ⱼ, ω, b::Vector{TT}, c::Vector{TT}, q::Vector{DT}, vᵢ::Vector{DT}, vⱼ::Vector{DT}, B::Matrix{DT}) where {DT,TT}
    @assert length(b) == length(c)

    local result = zero(Double64)

    for i in eachindex(b)
        for j in eachindex(b)
            γ(c[i], c[j], q)
            γ̇ᵢ(c[i], c[j], vᵢ)
            γ̇ⱼ(c[i], c[j], vⱼ)
            ω(t, q, B)

            result += b[i] * b[j] * vector_matrix_vector_product(vⱼ, B, vᵢ)
        end
    end

    return result
end

function surface_integral(t, x::AbstractMatrix{DT}, ω, nx, ny) where {DT}
    local b = [0.5, 0.5]
    local c = [0.0, 1.0]

    local q  = zeros(DT, size(x,1))
    local vᵢ = zeros(DT, size(x,1))
    local vⱼ = zeros(DT, size(x,1))
    local B  = zeros(DT, size(x,1), size(x,1))
    local I  = zero(Double64)

    integrate_trapezoidal = (γ, γ̇ᵢ, γ̇ⱼ) -> integrate(t, γ, γ̇ᵢ, γ̇ⱼ, ω, b, c, q, vᵢ, vⱼ, B)

    for j in 1:ny-1
        for i in 1:nx-1
            γ  = (λ, μ, y) -> interpolate_trajectory(x, i, j, λ, μ, y, nx, ny)
            γ̇ᵢ = (λ, μ, y) -> interpolate_derivative_i(x, i, j, λ, μ, y, nx, ny)
            γ̇ⱼ = (λ, μ, y) -> interpolate_derivative_j(x, i, j, λ, μ, y, nx, ny)
            I += integrate_trapezoidal(γ, γ̇ᵢ, γ̇ⱼ)
        end
    end

    return I
end


function integrate_canonical(γ̇ᵢ, γ̇ⱼ, Θ̇ᵢ, Θ̇ⱼ, b::Vector{TT}, c::Vector{TT}, vᵢ::Vector{DT}, vⱼ::Vector{DT}) where {DT,TT}
    @assert length(b) == length(c)

    local result = zero(Double64)

    for i in eachindex(b)
        for j in eachindex(b)
            Θ̇ᵢ(c[i], c[j], vᵢ)
            γ̇ⱼ(c[i], c[j], vⱼ)
            result += b[i] * b[j] * vᵢ ⋅ vⱼ

            γ̇ᵢ(c[i], c[j], vᵢ)
            Θ̇ⱼ(c[i], c[j], vⱼ)
            result -= b[i] * b[j] * vᵢ ⋅ vⱼ
        end
    end

    return result
end

function surface_integral_canonical(q::AbstractMatrix{DT}, p::AbstractMatrix{DT}, nx, ny) where {DT}
    local b = [0.5, 0.5]
    local c = [0.0, 1.0]

    local vᵢ = zeros(DT, size(q,1))
    local vⱼ = zeros(DT, size(q,1))
    local I  = zero(Double64)

    integrate_trapezoidal = (γ̇ᵢ, γ̇ⱼ, Θ̇ᵢ, Θ̇ⱼ) -> integrate_canonical(γ̇ᵢ, γ̇ⱼ, Θ̇ᵢ, Θ̇ⱼ, b, c, vᵢ, vⱼ)

    for j in 1:ny-1
        for i in 1:nx-1
            γ̇ᵢ = (λ, μ, y) -> interpolate_derivative_i(q, i, j, λ, μ, y, nx, ny)
            γ̇ⱼ = (λ, μ, y) -> interpolate_derivative_j(q, i, j, λ, μ, y, nx, ny)
            Θ̇ᵢ = (λ, μ, y) -> interpolate_derivative_i(p, i, j, λ, μ, y, nx, ny)
            Θ̇ⱼ = (λ, μ, y) -> interpolate_derivative_j(p, i, j, λ, μ, y, nx, ny)
            I += integrate_trapezoidal(γ̇ᵢ, γ̇ⱼ, Θ̇ᵢ, Θ̇ⱼ)
        end
    end

    return I
end


function evaluate_poincare_invariant(pinv::PoincareInvariant2ndTrapezoidal{DT}, sol::Solution) where {DT}
    local verbosity = get_config(:verbosity)

    for i in axes(sol.q,1)
        verbosity > 1 ? println("      it = ", i) : nothing
        @views pinv.I[i]  = surface_integral(sol.t[i], hcat(sol.q[i,:]...), pinv.ω, pinv.nx, pinv.ny)
        pinv.ΔI[i] = abs(pinv.I[0]) < sqrt(eps()) ? pinv.I[i] : (pinv.I[i] .- pinv.I[0]) ./ pinv.I[0]
        verbosity > 1 ? println("           I_q = ", pinv.I[i], ",   ε_q = ", pinv.ΔI[i]) : nothing

        if hasproperty(sol, :p)
            @views pinv.J[i] = surface_integral_canonical(hcat(sol.q[i,:]...), hcat(sol.p[i,:]...), pinv.nx, pinv.ny)
            pinv.ΔJ[i] = abs(pinv.J[0]) < sqrt(eps()) ? pinv.J[i] : (pinv.J[i] .- pinv.J[0]) ./ pinv.J[0]
            verbosity > 1 ? println("           I_p = ", pinv.J[i], ",   ε_p = ", pinv.ΔJ[i]) : nothing
        end
    end

    return (DT.(pinv.I), DT.(pinv.J), DT.(pinv.ΔI), DT.(pinv.ΔJ))
end


function write_to_hdf5(pinv::PoincareInvariant2ndTrapezoidal, sol::Solution, output_file::String)
    # h5open(output_file, isfile(output_file) ? "r+" : "w") do h5
    h5open(output_file, "w") do h5

        write(h5, "t", sol.t)
        write(h5, "I", pinv.I)

        isdefined(sol, :p) ? write(h5, "J", pinv.J) : nothing

    end
end

=#

end  # FiniteDifferences
