
import FastTransforms: paduapoints, plan_chebyshevtransform, plan_ichebyshevtransform
import OrthogonalPolynomialsQuasi

const OPQ = OrthogonalPolynomialsQuasi


struct PoincareInvariant2ndOPQ{DT,TT,ET,ΩT,CT} <: AbstractPoincareInvariant2nd{DT}
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
    y::Array{DT,3}
    z̃::Array{DT,3}
    z::Array{DT,3}
    Ω::Matrix{DT}
    cheb::CT
end

function PoincareInvariant2ndOPQ(f_equ::Function, f_surface::Function, ω::ΩT, Δt::TT, d::Int, nx::Int, ny::Int, ntime::Int, nsave::Int=1, DT=Float64) where {TT,ΩT}

    nx = 2nx
    ny = 2ny

    if get_config(:verbosity) > 1
        println()
        println("Second Euler-Poincaré Integral Invariant (OrthogonalPolynomials)")
        println("================================================================")
        println()
        println(" nx    = ", nx)
        println(" ny    = ", ny)
        println(" ntime = ", ntime)
        println(" nsave = ", nsave)
        println(" Δt    = ", Δt)
        println()
    end

    # compute initial conditions
    xGrid = OPQ.grid(OPQ.ChebyshevT()[:,1:nx])
    yGrid = OPQ.grid(OPQ.ChebyshevT()[:,1:ny])

    q₀ = reshape([f_surface((xGrid[i] + 1) / 2, (yGrid[j] + 1) / 2) for i in 1:nx, j in 1:ny], nx*ny)

    equ = f_equ(q₀)


    # compute initial conditions
    # _padua_length(N) = Int(cld(-3+sqrt(1+8N),2))

    # c = paduapoints(DT, _padua_length(nx*ny))

    # println(length(c))

    # q₀ = zeros(DT, (d, length(c)))

    # for i in axes(c,1)
    #     q₀[:,i] .= f_surface(c[i][1], c[i][2])
    # end

    # equ = f_equ(q₀)
    

    # create arrays for results
    nt = div(ntime, nsave)

    I  = OffsetArray(zeros(Double64, nt+1), 0:nt)
    J  = OffsetArray(zeros(Double64, nt+1), 0:nt)
    ΔI = OffsetArray(zeros(Double64, nt+1), 0:nt)
    ΔJ = OffsetArray(zeros(Double64, nt+1), 0:nt)

    y  = zeros(DT, d, nx,   ny  )
    z̃  = zeros(DT, d, nx,   ny-1)
    z  = zeros(DT, d, nx-1, ny-1)
    Ω  = zeros(DT, d, d)

    Di = zeros(DT, d, nx-1, ny  )
    Dj = zeros(DT, d, nx,   ny-1)
    D1 = zeros(DT, d, nx-1, ny-1)
    D2 = zeros(DT, d, nx-1, ny-1)

    Tx = OPQ.ChebyshevT()[:,1:nx]
    Ty = OPQ.ChebyshevT()[:,1:ny]
    Ux = OPQ.ChebyshevU()[:,1:nx-1]
    Uy = OPQ.ChebyshevU()[:,1:ny-1]

    Dx = OPQ.Derivative(axes(Tx,1)) * Tx
    Dy = OPQ.Derivative(axes(Ty,1)) * Ty

    PTx = plan_chebyshevtransform(zeros(nx), Val(1))
    PTy = plan_chebyshevtransform(zeros(ny), Val(1))

    iPUx = plan_ichebyshevtransform(zeros(nx-1), Val(2))
    iPUy = plan_ichebyshevtransform(zeros(ny-1), Val(2))

    PTUx = Ux \ Tx
    PTUy = Uy \ Ty

    cheb = (Di = Di, Dj = Dj, D1 = D1, D2 = D2,
            Tx = Tx, Ty = Ty, Ux = Ux, Uy = Uy,
            Dx = Dx, Dy = Dy,
            PTx  = PTx,  PTy  = PTy,
            iPUx = iPUx, iPUy = iPUy,
            PTUx = PTUx, PTUy = PTUy)

    PoincareInvariant2ndOPQ{DT,TT,typeof(equ),ΩT,typeof(cheb)}(equ, ω, Δt, nx, ny, ntime, nsave, nt, I, J, ΔI, ΔJ, y, z̃, z, Ω, cheb)
end



# function compute_position(x::AbstractArray{DT}, nx, ny, Ux, Uy, Px, Py, PUx, PUy, UxGrid, UyGrid, PUxi, PUyi) where {DT}
#     local y = copy(reshape(x, (size(x,1), nx, ny)))
#     local z̃ = zeros(DT, (size(x,1), nx,   ny-1))
#     local z = zeros(DT, (size(x,1), nx-1, ny-1))

#     for d in axes(y,1)
#         for i in axes(y,2)
#             @views y[d,i,:] .= Py \ y[d,i,:]
#             # @views y[d,i,:] .= Py * y[d,i,:]
#         end
#         for j in axes(y,3)
#             @views y[d,:,j] .= Px \ y[d,:,j]
#             # @views y[d,:,j] .= Px * y[d,:,j]
#         end
#     end

#     for d in axes(z̃,1)
#         for i in axes(z̃,2)
#             @views z̃[d,i,:] .= PUy * y[d,i,:]
#         end
#     end
#     for d in axes(z,1)
#         for j in axes(z,3)
#             @views z[d,:,j] .= PUx * z̃[d,:,j]
#         end
#     end

#     for d in axes(z,1)
#         for i in axes(z,2)
#             @views z[d,i,:] .= (Uy * z[d,i,:])[UyGrid]
#             # @views z[d,i,:] .= PUyi * z[d,i,:]
#         end
#         for j in axes(z,3)
#             @views z[d,:,j] .= (Ux * z[d,:,j])[UxGrid]
#             # @views z[d,:,j] .= PUxi * z[d,:,j]
#         end
#     end
    
#     return y, z
# end


# function compute_derivative_i!(Di, D1, y, Dx, Ux, Uy, PUy, UxGrid, UyGrid, PUxi, PUyi)
#     for d in axes(Di,1)
#         for j in axes(Di,3)
#             # @views Di[d,:,j] .= Ux \ ( Dx * y[d,:,j] )
#             @views Di[d,:,j] .= ( Dx * y[d,:,j] ).args[2].args[1]
            
#         end
#     end

#     for d in axes(D1,1)
#         for i in axes(D1,2)
#             @views D1[d,i,:] .= PUy * Di[d,i,:]
#         end
#     end

#     for d in axes(D1,1)
#         for i in axes(D1,2)
#             @views D1[d,i,:] .= (Uy * D1[d,i,:])[UyGrid]
#             # @views D1[d,i,:] .= PUyi * D1[d,i,:]
#         end
#         for j in axes(D1,3)
#             @views D1[d,:,j] .= (Ux * D1[d,:,j])[UxGrid]
#             # @views D1[d,:,j] .= PUxi * D1[d,:,j]
#         end
#     end
# end


# function compute_derivative_j!(Dj, D2, y, Dy, Ux, Uy, PUx, UxGrid, UyGrid, PUxi, PUyi)
#     for d in axes(Dj,1)
#         for i in axes(Dj,2)
# #            @views Dj[d,i,:] .= Uy \ ( Dy * y[d,i,:] )
#             @views Dj[d,i,:] .= ( Dy * y[d,i,:] ).args[2].args[1]
#         end
#     end

#     for d in axes(D2,1)
#         for j in axes(D2,3)
#             @views D2[d,:,j] .= PUx * Dj[d,:,j]
#         end
#     end

#     for d in axes(D2,1)
#         for i in axes(D2,2)
#             @views D2[d,i,:] .= (Uy * D2[d,i,:])[UyGrid]
#             # @views D2[d,i,:] .= PUyi * D2[d,i,:]
#         end
#         for j in axes(D2,3)
#             @views D2[d,:,j] .= (Ux * D2[d,:,j])[UxGrid]
#             # @views D2[d,:,j] .= PUxi * D2[d,:,j]
#         end
#     end
# end



# function compute_position(x::AbstractArray{DT}, nx, ny, PTx, PTy, PTUx, PTUy, iPUx, iPUy) where {DT}
#     local y = copy(reshape(x, (size(x,1), nx, ny)))
#     local z̃ = zeros(DT, (size(x,1), nx,   ny-1))
#     local z = zeros(DT, (size(x,1), nx-1, ny-1))

#     # Py = iPUy * PTUy * PTy
#     # Px = iPUx * PTUx * PTx

#     for d in axes(z̃,1)
#         for i in axes(z̃,2)
#             # @views z̃[d,i,:] .= Py * y[d,i,:]
#             @views z̃[d,i,:] .= iPUy * ( PTUy * ( PTy * y[d,i,:] ) )
#             # @views z̃[d,i,:] .= iPUy * ( PTUy * ( PTy \ y[d,i,:] ) )
#         end
#     end

#     for d in axes(z,1)
#         for j in axes(z,3)
#             # @views z[d,:,j] .= Px * z̃[d,:,j]
#             @views z[d,:,j] .= iPUx * ( PTUx * ( PTx * z̃[d,:,j] ) )
#             # @views z[d,:,j] .= iPUx * ( PTUx * ( PTx \ z̃[d,:,j] ) )
#         end
#     end
    
#     return z
# end


# function compute_derivative_i!(Di, D1, x, nx, ny, Dx, PTx, PTy, PTUy, iPUx, iPUy)
#     local y = copy(reshape(x, (size(x,1), nx, ny)))

#     # Py = iPUy * PTUy * PTy

#     for d in axes(Di,1)
#         for j in axes(Di,3)
#             @views Di[d,:,j] .= iPUx * ( Dx * (PTx * y[d,:,j]) ).args[2].args[1]
#             # @views Di[d,:,j] .= iPUx * ( Dx * (PTx \ y[d,:,j]) ).args[2].args[1]
#         end
#     end

#     for d in axes(D1,1)
#         for i in axes(D1,2)
#             # @views D1[d,i,:] .= Py * Di[d,i,:]
#             @views D1[d,i,:] .= iPUy * ( PTUy * ( PTy * Di[d,i,:] ) )
#             # @views D1[d,i,:] .= iPUy * ( PTUy * ( PTy \ Di[d,i,:] ) )
#         end
#     end
# end


# function compute_derivative_j!(Dj, D2, x, nx, ny, Dy, PTx, PTy, PTUx, iPUx, iPUy)
#     local y = copy(reshape(x, (size(x,1), nx, ny)))

#     # Px = iPUx * PUx * PTx

#     for d in axes(Dj,1)
#         for i in axes(Dj,2)
#             @views Dj[d,i,:] .= iPUy * ( Dy * (PTy * y[d,i,:]) ).args[2].args[1]
#             # @views Dj[d,i,:] .= iPUy * ( Dy * (PTy \ y[d,i,:]) ).args[2].args[1]
#         end
#     end

#     for d in axes(D2,1)
#         for j in axes(D2,3)
#             # @views D2[d,:,j] .= Px * Dj[d,:,j]
#             @views D2[d,:,j] .= iPUx * ( PTUx * ( PTx * Dj[d,:,j] ) )
#             # @views D2[d,:,j] .= iPUx * ( PTUx * ( PTx \ Dj[d,:,j] ) )
#         end
#     end
# end


# function _integrate_opq(t, ω, z, nx, ny, D1, D2, Ω::Matrix)
#     local result::Double64 = 0

#     for j in axes(z,3)
#         for i in axes(z,2)
#             ω(t, z[:,i,j], Ω)
#             @views result += D2[:,i,j]' * Ω * D1[:,i,j]
#         end
#     end
    
#     return 4 * result / nx / ny
# end



function _integrate_opq(t, x, pinv)
    Di, Dj, D1, D2, Tx, Ty, Ux, Uy, Dx, Dy, PTx, PTy, iPUx, iPUy, PTUx, PTUy = pinv.cheb

    pinv.y .= reshape(x, (size(x,1), pinv.nx, pinv.ny))

    for d in axes(pinv.y,1)
        for i in axes(pinv.y,2)
            @views pinv.y[d,i,:] .= PTy * pinv.y[d,i,:]
        end
        for j in axes(pinv.y,3)
            @views pinv.y[d,:,j] .= PTx * pinv.y[d,:,j]
        end
    end

    for d in axes(pinv.z̃,1)
        for i in axes(pinv.z̃,2)
            @views pinv.z̃[d,i,:] .= PTUy * pinv.y[d,i,:]
        end
    end
    for d in axes(pinv.z,1)
        for j in axes(pinv.z,3)
            @views pinv.z[d,:,j] .= PTUx * pinv.z̃[d,:,j]
        end
    end

    for d in axes(pinv.z,1)
        for i in axes(pinv.z,2)
            @views pinv.z[d,i,:] .= iPUy * pinv.z[d,i,:]
        end
        for j in axes(pinv.z,3)
            @views pinv.z[d,:,j] .= iPUx * pinv.z[d,:,j]
        end
    end


    for d in axes(Di,1)
        for j in axes(Di,3)
            @views Di[d,:,j] .= ( Dx * pinv.y[d,:,j] ).args[2].args[1]
        end
    end
    for d in axes(D1,1)
        for i in axes(D1,2)
            @views D1[d,i,:] .= iPUy * (PTUy * Di[d,i,:])
        end
        for j in axes(D1,3)
            @views D1[d,:,j] .= iPUx * D1[d,:,j]
        end
    end

    for d in axes(Dj,1)
        for i in axes(Dj,2)
            @views Dj[d,i,:] .= ( Dy * pinv.y[d,i,:] ).args[2].args[1]
        end
    end
    for d in axes(D2,1)
        for j in axes(D2,3)
            @views D2[d,:,j] .= iPUx * (PTUx * Dj[d,:,j])
        end
        for i in axes(D2,2)
            @views D2[d,i,:] .= iPUy * D2[d,i,:]
        end
    end


    local result::Double64 = 0

    for j in axes(pinv.z,3)
        for i in axes(pinv.z,2)
            @views pinv.ω(t, pinv.z[:,i,j], pinv.Ω)
            @views result += D2[:,i,j]' * pinv.Ω * D1[:,i,j]
        end
    end
    
    return 4 * result / pinv.nx / pinv.ny
end


function surface_integral(pinv, t, x::AbstractMatrix{DT}) where {DT}
    # local B = zeros(DT, size(x,1), size(x,1))
    # local I = zero(Double64)


    # y, z = compute_position(x, pinv.nx, pinv.ny, Ux, Uy, Px, Py, PUx, PUy, UxGrid, UyGrid)
    # y, z = compute_position(x, pinv.nx, pinv.ny, Ux, Uy, Px, Py, PUx, PUy, UxGrid, UyGrid, UxiPlan, UyiPlan)

    # z = compute_position(x, pinv.nx, pinv.ny, PTx, PTy, PTUx, PTUy, iPUx, iPUy)

    # compute_derivative_i!(Di, D1, x, pinv.nx, pinv.ny, Dx, PTx, PTy, PTUy, iPUx, iPUy)
    # compute_derivative_j!(Dj, D2, x, pinv.nx, pinv.ny, Dy, PTx, PTy, PTUx, iPUx, iPUy)

    # return _integrate_opq(t, pinv.ω, z, pinv.nx, pinv.ny, D1, D2, B)

    return _integrate_opq(t, x, pinv)
end


# function integrate_canonical(γ̇ᵢ, γ̇ⱼ, Θ̇ᵢ, Θ̇ⱼ, b::Vector{TT}, c::Vector{TT}, vᵢ::Vector{DT}, vⱼ::Vector{DT}) where {DT,TT}
#     @assert length(b) == length(c)

#     local result = zero(Double64)

#     for i in eachindex(b)
#         for j in eachindex(b)
#             Θ̇ᵢ(c[i], c[j], vᵢ)
#             γ̇ⱼ(c[i], c[j], vⱼ)
#             result += b[i] * b[j] * dot(vᵢ,vⱼ)

#             γ̇ᵢ(c[i], c[j], vᵢ)
#             Θ̇ⱼ(c[i], c[j], vⱼ)
#             result -= b[i] * b[j] * dot(vᵢ,vⱼ)
#         end
#     end

#     return result
# end


# function surface_integral_canonical(q::AbstractMatrix{DT}, p::AbstractMatrix{DT}, nx, ny) where {DT}
#     local b = [0.5, 0.5]
#     local c = [0.0, 1.0]

#     local vᵢ = zeros(DT, size(q,1))
#     local vⱼ = zeros(DT, size(q,1))
#     local I  = zero(Double64)

#     integrate_trapezoidal = (γ, γ̇ᵢ, γ̇ⱼ, Θ̇ᵢ, Θ̇ⱼ) -> integrate_canonical(γ̇ᵢ, γ̇ⱼ, Θ̇ᵢ, Θ̇ⱼ, b, c, vᵢ, vⱼ)

#     for j in 1:ny-1
#         for i in 1:nx-1
#             γ̇ᵢ = (λ, μ, y) -> interpolate_derivative_i(q, i, j, λ, μ, y, nx, ny)
#             γ̇ⱼ = (λ, μ, y) -> interpolate_derivative_j(q, i, j, λ, μ, y, nx, ny)
#             Θ̇ᵢ = (λ, μ, y) -> interpolate_derivative_i(p, i, j, λ, μ, y, nx, ny)
#             Θ̇ⱼ = (λ, μ, y) -> interpolate_derivative_j(p, i, j, λ, μ, y, nx, ny)
#             I += integrate_trapezoidal(γ, γ̇ᵢ, γ̇ⱼ, Θ̇ᵢ, Θ̇ⱼ)
#         end
#     end

#     return I
# end


function evaluate_poincare_invariant(pinv::PoincareInvariant2ndOPQ{DT}, sol::Solution) where {DT}
    local verbosity = get_config(:verbosity)

    for i in axes(sol.q,2)
        verbosity > 1 ? println("      it = ", i) : nothing
        @views pinv.I[i]  = surface_integral(pinv, sol.t[i], hcat(sol.q[i,:]...))
        pinv.ΔI[i] = abs(pinv.I[0]) < sqrt(eps()) ? pinv.I[i] : (pinv.I[i] .- pinv.I[0]) ./ pinv.I[0]
        verbosity > 1 ? println("           I_q = ", pinv.I[i], ",   ε_q = ", pinv.ΔI[i]) : nothing

        # if hasproperty(sol, :p)
        #     pinv.J[i] = surface_integral_canonical(pinv, hcat(sol.q[i,:]...), hcat(sol.p[i,:]...))
        #     pinv.ΔJ[i] = abs(pinv.J[0]) < sqrt(eps()) ? pinv.J[i] : (pinv.J[i] .- pinv.J[0]) ./ pinv.J[0]
        #     verbosity > 1 ? println("           I_p = ", pinv.J[i], ",   ε_p = ", pinv.ΔJ[i]) : nothing
        # end
    end

    return (DT.(pinv.I), DT.(pinv.J), DT.(pinv.ΔI), DT.(pinv.ΔJ))
end


function write_to_hdf5(pinv::PoincareInvariant2ndOPQ, sol::Solution, output_file::String)
    # h5open(output_file, isfile(output_file) ? "r+" : "w") do h5
    h5open(output_file, "w") do h5

        write(h5, "t", sol.t)
        write(h5, "I", pinv.I)

        isdefined(sol, :p) ? write(h5, "J", pinv.J) : nothing

    end
end
