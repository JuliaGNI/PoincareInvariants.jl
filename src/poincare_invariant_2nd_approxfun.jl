
import ApproxFun
import ApproxFun: Fun, ⊗, (..)

## 2nd Poincaré Invariant for general two-form

struct PoincareInvariant2ndApproxFun{DT,ND,NC,NV,ET,ΩT,ϑT} <: AbstractPoincareInvariant2ndApproxFun{DT,ND,NC,NV}
    equ::ET
    ω::ΩT
    D²ϑ::ϑT
    Δt::DT
    nx::Int
    ny::Int
    ntime::Int
    nsave::Int
    nt::Int
    I::OffsetArray{Double64,1,Vector{Double64}}
    J::OffsetArray{Double64,1,Vector{Double64}}
    K::OffsetArray{Double64,1,Vector{Double64}}
    L::OffsetArray{Double64,1,Vector{Double64}}
    ΔI::OffsetArray{Double64,1,Vector{Double64}}
    ΔJ::OffsetArray{Double64,1,Vector{Double64}}
end

function PoincareInvariant2ndApproxFun(f_equ::Function, f_surface::Function, ω::ΩT, D²ϑ::ϑT, Δt::TT, nd::Int, nx::Int, ny::Int, ntime::Int, nsave::Int=1, DT=Float64) where {TT,ΩT,ϑT}

    # compute Chebyshev points
    c = ApproxFun.points(ApproxFun.Chebyshev(0..1)^2, nx*ny)

    # compute initial conditions
    q₀ = [f_surface(c[i][1], c[i][2]) for i in eachindex(c)]

    # initialise euation
    equ = f_equ(q₀)

    # create arrays for results
    nt = div(ntime, nsave)

    I  = OffsetArray(zeros(Double64, nt+1), 0:nt)
    J  = OffsetArray(zeros(Double64, nt+1), 0:nt)
    K  = OffsetArray(zeros(Double64, nt+1), 0:nt)
    L  = OffsetArray(zeros(Double64, nt+1), 0:nt)
    ΔI = OffsetArray(zeros(Double64, nt+1), 0:nt)
    ΔJ = OffsetArray(zeros(Double64, nt+1), 0:nt)

    # get size of coefficient and value vectors
    SC = ApproxFun.Chebyshev(0..1)^2
    SU = ApproxFun.Ultraspherical(1, 0..1)^2
    Dx = ApproxFun.Derivative(SC, [1,0])
    Dy = ApproxFun.Derivative(SC, [0,1])

    fqc = Fun(SC, ApproxFun.transform(SC, hcat(q₀...)[1,:]))
    fqu = Fun(Dx * fqc, SU)
    nc = ApproxFun.ncoefficients(fqu)
    nv = length(ApproxFun.values(fqu))

    if get_config(:verbosity) > 1
        println()
        println("Second Euler-Poincaré Integral Invariant (ApproxFun)")
        println("====================================================")
        println()
        println(" nx    = ", nx)
        println(" ny    = ", ny)
        println(" np    = ", length(c))
        println(" nc(CC)= ", ApproxFun.ncoefficients(fqc))
        println(" nv(CC)= ", length(ApproxFun.values(fqc)))
        println(" nc(US)= ", nc)
        println(" nv(US)= ", nv)
        println(" ntime = ", ntime)
        println(" nsave = ", nsave)
        println(" Δt    = ", Δt)
        println()
    end

    # initialise Poincare invariant
    PoincareInvariant2ndApproxFun{DT,nd,nc,nv,typeof(equ),ΩT,ϑT}(equ, ω, D²ϑ, DT(Δt), nx, ny, ntime, nsave, nt, I, J, K, L, ΔI, ΔJ)
end


function PoincareInvariant2ndApproxFun(f_equ::Function, f_surface::Function, ω::ΩT, Δt::TT, nd::Int, nx::Int, ny::Int, ntime::Int, nsave::Int=1, DT=Float64) where {TT,ΩT}
    PoincareInvariant2ndApproxFun(f_equ, f_surface, ω, (), Δt, nd, nx, ny, ntime, nsave, DT)
end


function evaluate_poincare_invariant(pinv::PoincareInvariant2ndApproxFun{DT}, sol::Solution) where {DT}

    local verbosity = get_config(:verbosity)

    # local q = permutedims(sol.q, [3,1,2])
    # if isdefined(sol, :p)
    #     local p = permutedims(sol.p, [3,1,2])
    # end
    # if isdefined(sol, :λ)
    #     local λ = permutedims(sol.λ, [3,1,2])
    # end

    verbosity ≤ 1 ? prog = Progress(size(sol.q,1), 5) : nothing

    for i in axes(sol.q,1)
        verbosity > 1 ? println("      it = ", i) : nothing
        pinv.I[i] = compute_noncanonical_invariant(pinv, sol.t[i], hcat(sol.q[i,:]...))
        pinv.ΔI[i] = abs(pinv.I[0]) < sqrt(eps()) ? pinv.I[i] : (pinv.I[i] .- pinv.I[0]) ./ pinv.I[0]
        verbosity > 1 ? println("           I_q = ", pinv.I[i], ",   ε_q = ", pinv.ΔI[i]) : nothing
    end

    if isdefined(sol, :p)
        for i in axes(sol.q,1)
            pinv.J[i] = compute_canonical_invariant(pinv, hcat(sol.q[i,:]...), hcat(sol.p[i,:]...))
            pinv.ΔJ[i] = abs(pinv.J[0]) < sqrt(eps()) ? pinv.J[i] : (pinv.J[i] .- pinv.J[0]) ./ pinv.J[0]
            verbosity > 1 ? println("           I_p = ", pinv.J[i], ",   ε_p = ", pinv.ΔJ[i]) : nothing
        end
    end

    return (DT.(pinv.I), DT.(pinv.J), DT.(pinv.ΔI), DT.(pinv.ΔJ))
end


function evaluate_poincare_invariant_correction(pinv::PoincareInvariant2ndApproxFun{DT}, sol::Solution) where {DT}

    local verbosity = get_config(:verbosity)

    if isdefined(sol, :λ)
        for i in axes(sol.q,1)
            pinv.K[i] = compute_noncanonical_correction(pinv, sol.t[i], hcat(sol.q[i,:]...), hcat(sol.λ[i,:]...))
            pinv.L[i] = pinv.I[i] - pinv.Δt^2 * pinv.K[i]
            verbosity > 1 ? println("           I_λ = ", pinv.L[i], ",   ε_λ = ", (pinv.L[i]-pinv.L[0])/pinv.L[0]) : nothing
            verbosity > 1 ? println("           K_λ = ", pinv.K[i]) : nothing
        end

        verbosity ≤ 1 ? next!(prog) : nothing
    end

    return (DT.(pinv.K), DT.(pinv.L))
end


function write_to_hdf5(pinv::PoincareInvariant2ndApproxFun{DT}, sol::Solution, output_file::String) where {DT}
    # h5open(output_file, isfile(output_file) ? "r+" : "w") do h5
    h5open(output_file, "w") do h5

        write(h5, "t", sol.t)
        write(h5, "I", pinv.I)

        isdefined(sol, :p) ? write(h5, "J", pinv.J) : nothing
        isdefined(sol, :λ) ? write(h5, "K", pinv.L) : nothing
        isdefined(sol, :λ) ? write(h5, "L", pinv.L) : nothing
    end
end


## 2nd Poincaré Invariant for canonical two-form

struct PoincareInvariant2ndApproxFunCanonical{DT,ND,NC,NV,ET} <: AbstractPoincareInvariant2ndApproxFun{DT,ND,NC,NV}
    equ::ET
    Δt::DT
    nx::Int
    ny::Int
    ntime::Int
    nsave::Int
    nt::Int
    I::OffsetArray{Double64,1,Vector{Double64}}
    ΔI::OffsetArray{Double64,1,Vector{Double64}}
end

function PoincareInvariant2ndApproxFunCanonical(f_equ::Function, f_surface_q::Function, f_surface_p::Function, Δt::TT, nd::Int, nx::Int, ny::Int, ntime::Int, nsave::Int=1, DT=Float64) where {TT}

    if get_config(:verbosity) > 1
        println()
        println("Second Canonical Euler-Poincaré Integral Invariant (ApproxFun)")
        println("==============================================================")
        println()
        println(" nx    = ", nx)
        println(" ny    = ", ny)
        println(" ntime = ", ntime)
        println(" nsave = ", nsave)
        println(" Δt    = ", Δt)
        println()
    end

    # compute Chebyshev points
    c = ApproxFun.points(ApproxFun.Chebyshev(0..1)^2, nx*ny)

    # compute initial conditions
    q₀ = [f_surface_q(c[i][1], c[i][2]) for i in eachindex(c)]
    p₀ = [f_surface_p(c[i][1], c[i][2]) for i in eachindex(c)]

    # initialise euation
    equ = f_equ(q₀, p₀)

    # create arrays for results
    nt = div(ntime, nsave)

    I  = OffsetArray(zeros(Double64, nt+1), 0:nt)
    ΔI = OffsetArray(zeros(Double64, nt+1), 0:nt)

    # get size of coefficient and value vectors
    SC = ApproxFun.Chebyshev(0..1)^2
    SU = ApproxFun.Ultraspherical(1, 0..1)^2
    Dx = ApproxFun.Derivative(SC, [1,0])
    Dy = ApproxFun.Derivative(SC, [0,1])

    fq = Fun(Dx * Fun(SC, ApproxFun.transform(SC, hcat(q₀...)[1,:])), SU)
    nc = ApproxFun.ncoefficients(fq)
    nv = length(ApproxFun.values(fq))

    # initialise Poincare invariant
    PoincareInvariant2ndApproxFunCanonical{DT,nd,nc,nv,typeof(equ)}(equ, DT(Δt), nx, ny, ntime, nsave, nt, I, ΔI)
end


function evaluate_poincare_invariant(pinv::PoincareInvariant2ndApproxFunCanonical{DT}, sol::Solution) where {DT}
    local verbosity = get_config(:verbosity)

    verbosity ≤ 1 ? prog = Progress(size(sol.q,1), 5) : nothing

    for i in axes(sol.q,1)
        verbosity > 1 ? println("      it = ", i) : nothing
        pinv.I[i] = compute_canonical_invariant(pinv, hcat(sol.q[i,:]...), hcat(sol.p[i,:]...))
        verbosity > 1 ? println("           I_p = ", pinv.I[i], ",   ε_p = ", (pinv.I[i]-pinv.I[0])/pinv.I[0]) : nothing

        verbosity ≤ 1 ? next!(prog) : nothing
    end

    return pinv.I
end


function write_to_hdf5(pinv::PoincareInvariant2ndApproxFunCanonical{DT}, sol::Solution, output_file::String) where {DT}
    # h5open(output_file, isfile(output_file) ? "r+" : "w") do h5
    h5open(output_file, "w") do h5
        write(h5, "t", sol.t)
        write(h5, "I", pinv.I)
    end
end


## Common functionality

@generated function compute_integral(Is)
    SU = ApproxFun.Ultraspherical(1, 0..1)^2
    Q  = ApproxFun.DefiniteIntegral(SU)

    quote
        return Number( $Q * Fun($SU, ApproxFun.transform($SU, Is)) )
    end
end

@generated function compute_derivative1(q)
    Dx = ApproxFun.Derivative(ApproxFun.Chebyshev(0..1)^2, [1,0])
    CU = ApproxFun.Conversion(ApproxFun.Ultraspherical(1, 0..1) ⊗ ApproxFun.Chebyshev(0..1), ApproxFun.Ultraspherical(1, 0..1)^2)

    quote
        return values($CU * ($Dx * q))
    end
end

@generated function compute_derivative2(q)
    Dy = ApproxFun.Derivative(ApproxFun.Chebyshev(0..1)^2, [0,1])
    CU = ApproxFun.Conversion(ApproxFun.Chebyshev(0..1) ⊗ ApproxFun.Ultraspherical(1, 0..1), ApproxFun.Ultraspherical(1, 0..1)^2)

    quote
        return values($CU * ($Dy * q))
    end
end

@generated function compute_canonical_invariant(pinv::AbstractPoincareInvariant2ndApproxFun{DT,ND,NC,NV}, q::AbstractArray{DT,2}, p::AbstractArray{DT,2}) where {DT,ND,NC,NV}
    SC = ApproxFun.Chebyshev(0..1)^2
    SU = ApproxFun.Ultraspherical(1, 0..1)^2

    Js = zeros(DT, NV)
    qσ = zeros(DT, ND, NV)
    qτ = zeros(DT, ND, NV)
    pσ = zeros(DT, ND, NV)
    pτ = zeros(DT, ND, NV)

    quote
        for j in axes(p,1)
            fq = Fun($SC, ApproxFun.transform($SC, q[j,:]))
            fp = Fun($SC, ApproxFun.transform($SC, p[j,:]))

            # compute derivatives of q and store values
            $qσ[j,:] .= compute_derivative1(fq)
            $qτ[j,:] .= compute_derivative2(fq)

            # compute derivatives of p and store values
            $pσ[j,:] .= compute_derivative1(fp)
            $pτ[j,:] .= compute_derivative2(fp)
        end

        # compute integrands of integral invariants at all points
        for k in 1:NV
            $Js[k] = $pσ[:,k] ⋅ $qτ[:,k] - $qσ[:,k] ⋅ $pτ[:,k]
        end

        # compute canonical integral invariant
        return compute_integral($Js)
    end
end

@generated function compute_noncanonical_invariant(pinv::PoincareInvariant2ndApproxFun{DT,ND,NC,NV}, t::DT, q::AbstractArray{DT,2}) where {DT,ND,NC,NV}
    SC = ApproxFun.Chebyshev(0..1)^2
    SU = ApproxFun.Ultraspherical(1, 0..1)^2
    CU = ApproxFun.Conversion(SC, SU)

    Is = zeros(DT, NV)
    qs = zeros(DT, ND, NV)
    qσ = zeros(DT, ND, NV)
    qτ = zeros(DT, ND, NV)
    Ω  = zeros(DT, ND, ND)

    quote
        # loop over dimensions of state vector
        for j in axes(q,1)
            # obtain Chebyshev coefficients
            fq = Fun($SC, ApproxFun.transform($SC, q[j,:]))

            # obtain values of q at the same points as the derivatives
            $qs[j,:] .= ApproxFun.itransform($SU, resize!(collect(($CU * fq).coefficients), NC))

            # compute derivatives of q and store values
            $qσ[j,:] .= compute_derivative1(fq)
            $qτ[j,:] .= compute_derivative2(fq)
        end

        # compute integrands of integral invariants at all points
        for k in 1:NV
            pinv.ω(t, $qs[:,k], $Ω)
            $Is[k] = vector_matrix_vector_product($qτ[:,k], $Ω, $qσ[:,k])
        end

        # compute noncanonical integral invariant
        return compute_integral($Is)
    end
end


@generated function compute_noncanonical_correction(pinv::PoincareInvariant2ndApproxFun{DT,ND,NC,NV}, t::DT, q::AbstractArray{DT,2}, λ::AbstractArray{DT,2}) where {DT,ND,NC,NV}
    SC = ApproxFun.Chebyshev(0..1)^2
    SU = ApproxFun.Ultraspherical(1, 0..1)^2
    CU = ApproxFun.Conversion(SC, SU)

    Ks = zeros(DT, NV)
    qs = zeros(DT, ND, NV)
    λs = zeros(DT, ND, NV)
    qσ = zeros(DT, ND, NV)
    qτ = zeros(DT, ND, NV)
    λσ = zeros(DT, ND, NV)
    λτ = zeros(DT, ND, NV)
    Ω  = zeros(DT, ND, ND)
    D²ϑ= zeros(DT, ND, ND)

    quote
        for j in axes(q,1)
            fq = Fun($SC, ApproxFun.transform($SC, q[j,:]))
            fλ = Fun($SC, ApproxFun.transform($SC, λ[j,:]))

            # obtain values of q at the same points as the derivatives
            $qs[j,:] .= ApproxFun.itransform($SU, resize!(($CU * fq).coefficients, NC))
            $λs[j,:] .= ApproxFun.itransform($SU, resize!(($CU * fλ).coefficients, NC))

            # compute derivatives of q and store values
            $qσ[j,:] .= compute_derivative1(fq)
            $qτ[j,:] .= compute_derivative2(fq)

            # compute derivatives of λ and store values
            $λσ[j,:] .= compute_derivative1(fλ)
            $λτ[j,:] .= compute_derivative2(fλ)
        end

        # compute integrands of integral invariants at all points
        for k in 1:NV
            pinv.ω(t, $qs[:,k], $Ω)
            $Ks[k] = vector_matrix_vector_product($λτ[:,k], $Ω, $λσ[:,k])

            for j in eachindex(pinv.D²ϑ)
                pinv.D²ϑ[j](t, $qs[:,k], $D²ϑ)
                $Ks[k] += 2 * $λs[j,k] * vector_matrix_vector_product($qτ[:,k], $D²ϑ, $λσ[:,k])
                $Ks[k] -= 2 * $λs[j,k] * vector_matrix_vector_product($qσ[:,k], $D²ϑ, $λτ[:,k])
            end
        end

        # compute correction to noncanonical integral invariant
        return compute_integral($Ks)
    end
end
