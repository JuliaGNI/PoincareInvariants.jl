
## 2nd Poincaré Invariant for general two-form

struct PoincareInvariant2nd{DT,ND,NC,NV,ET,ΩT,ϑT}
    equ::ET
    ω::ΩT
    D²ϑ::ϑT
    Δt::DT
    nx::Int
    ny::Int
    ntime::Int
    nsave::Int
    nt::Int
    I::OffsetArray{DT,1,Vector{DT}}
    J::OffsetArray{DT,1,Vector{DT}}
    K::OffsetArray{DT,1,Vector{DT}}
    L::OffsetArray{DT,1,Vector{DT}}
end

function PoincareInvariant2nd(f_equ::Function, f_surface::Function, ω::ΩT, D²ϑ::ϑT, Δt::TT, nd::Int, nx::Int, ny::Int, ntime::Int, nsave::Int=1, DT=Float64) where {TT,ΩT,ϑT}

    if get_config(:verbosity) > 1
        println()
        println("Second Euler-Poincaré Integral Invariant")
        println("========================================")
        println()
        println(" nx    = ", nx)
        println(" ny    = ", ny)
        println(" ntime = ", ntime)
        println(" nsave = ", nsave)
        println(" Δt    = ", Δt)
        println()
    end

    # compute Chebyshev points
    c = points(Chebyshev(0..1)^2, nx*ny)

    # compute initial conditions
    q₀ = zeros(DT, (nd, length(c)))

    for i in eachindex(c)
        q₀[:,i] .= f_surface(c[i][1], c[i][2])
    end

    # initialise euation
    equ = f_equ(q₀)

    # create arrays for results
    nt = div(ntime, nsave)

    I = OffsetArray(zeros(DT, nt+1), 0:nt)
    J = OffsetArray(zeros(DT, nt+1), 0:nt)
    K = OffsetArray(zeros(DT, nt+1), 0:nt)
    L = OffsetArray(zeros(DT, nt+1), 0:nt)

    # get size of coefficient and value vectors
    SC = Chebyshev(0..1)^2
    SU = Ultraspherical(1, 0..1)^2
    Dx = Derivative(SC, [1,0])
    Dy = Derivative(SC, [0,1])

    fq = Fun(Dx * Fun(SC, ApproxFun.transform(SC, q₀[1,:])), SU)
    nc = ncoefficients(fq)
    nv = length(values(fq))

    # initialise Poincare invariant
    PoincareInvariant2nd{DT,nd,nc,nv,typeof(equ),ΩT,ϑT}(equ, ω, D²ϑ, DT(Δt), nx, ny, ntime, nsave, nt, I, J, K, L)
end


function PoincareInvariant2nd(f_equ::Function, f_surface::Function, ω::ΩT, Δt::TT, nd::Int, nx::Int, ny::Int, ntime::Int, nsave::Int=1, DT=Float64) where {TT,ΩT}
    PoincareInvariant2nd(f_equ, f_surface, ω, (), Δt, nd, nx, ny, ntime, nsave, DT)
end


function evaluate_poincare_invariant(pinv::PoincareInvariant2nd{DT}, sol::Solution) where {DT}

    local verbosity = get_config(:verbosity)

    # local q = permutedims(sol.q, [3,1,2])
    # if isdefined(sol, :p)
    #     local p = permutedims(sol.p, [3,1,2])
    # end
    # if isdefined(sol, :λ)
    #     local λ = permutedims(sol.λ, [3,1,2])
    # end

    verbosity ≤ 1 ? prog = Progress(size(sol.q,2), 5) : nothing

    for i in axes(sol.q,2)
        verbosity > 1 ? println("      it = ", i-1) : nothing
        pinv.I[i] = compute_noncanonical_invariant(pinv, sol.t[i], sol.q[:,i,:])
        verbosity > 1 ? println("           I_q = ", pinv.I[i], ",   ε_q = ", (pinv.I[i]-pinv.I[0])/pinv.I[0]) : nothing

        if isdefined(sol, :p)
            pinv.J[i] = compute_canonical_invariant(pinv, sol.q[:,i,:], sol.p[:,i,:])
            verbosity > 1 ? println("           I_p = ", pinv.J[i], ",   ε_p = ", (pinv.J[i]-pinv.J[0])/pinv.J[0]) : nothing
        end

        if isdefined(sol, :λ)
            pinv.K[i] = compute_noncanonical_correction(pinv, sol.t[i], sol.q[:,i,:], sol.λ[:,i,:])
            pinv.L[i] = pinv.I[i] - pinv.Δt^2 * pinv.K[i]
            verbosity > 1 ? println("           I_λ = ", pinv.L[i], ",   ε_λ = ", (pinv.L[i]-pinv.L[0])/pinv.L[0]) : nothing
            verbosity > 1 ? println("           K_λ = ", pinv.K[i]) : nothing
        end

        verbosity ≤ 1 ? next!(prog) : nothing
    end

    return (pinv.I, pinv.J, pinv.K, pinv.L)
end


function CommonFunctions.write_to_hdf5(pinv::PoincareInvariant2nd{DT}, sol::Solution, output_file::String) where {DT}
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

struct PoincareInvariant2ndCanonical{DT,ND,NC,NV,ET}
    equ::ET
    Δt::DT
    nx::Int
    ny::Int
    ntime::Int
    nsave::Int
    nt::Int
    I::OffsetArray{DT,1,Vector{DT}}
end

function PoincareInvariant2ndCanonical(f_equ::Function, f_surface_q::Function, f_surface_p::Function, Δt::TT, nd::Int, nx::Int, ny::Int, ntime::Int, nsave::Int=1, DT=Float64) where {TT}

    if get_config(:verbosity) > 1
        println()
        println("Second Euler-Poincaré Integral Invariant")
        println("========================================")
        println()
        println(" nx    = ", nx)
        println(" ny    = ", ny)
        println(" ntime = ", ntime)
        println(" nsave = ", nsave)
        println(" Δt    = ", Δt)
        println()
    end

    # compute Chebyshev points
    c = points(Chebyshev(0..1)^2, nx*ny)

    # compute initial conditions
    q₀ = zeros(DT, (nd, length(c)))
    p₀ = zeros(DT, (nd, length(c)))

    for i in eachindex(c)
        q₀[:,i] .= f_surface_q(c[i][1], c[i][2])
        p₀[:,i] .= f_surface_p(c[i][1], c[i][2])
    end

    # initialise euation
    equ = f_equ(q₀, p₀)

    # create arrays for results
    nt = div(ntime, nsave)

    I = OffsetArray(zeros(DT, nt+1), 0:nt)

    # get size of coefficient and value vectors
    SC = Chebyshev(0..1)^2
    SU = Ultraspherical(1, 0..1)^2
    Dx = Derivative(SC, [1,0])
    Dy = Derivative(SC, [0,1])

    fq = Fun(Dx * Fun(SC, ApproxFun.transform(SC, q₀[1,:])), SU)
    nc = ncoefficients(fq)
    nv = length(values(fq))

    # initialise Poincare invariant
    PoincareInvariant2ndCanonical{DT,nd,nc,nv,typeof(equ)}(equ, DT(Δt), nx, ny, ntime, nsave, nt, I)
end


function evaluate_poincare_invariant(pinv::PoincareInvariant2ndCanonical{DT}, sol::Solution) where {DT}
    local verbosity = get_config(:verbosity)

    verbosity ≤ 1 ? prog = Progress(size(sol.q,2), 5) : nothing

    for i in axes(sol.q,2)
        verbosity > 1 ? println("      it = ", i) : nothing
        pinv.I[i] = compute_canonical_invariant(pinv, sol.q[:,i,:], sol.p[:,i,:])
        verbosity > 1 ? println("           I_p = ", pinv.I[i], ",   ε_p = ", (pinv.I[i]-pinv.I[0])/pinv.I[0]) : nothing

        verbosity ≤ 1 ? next!(prog) : nothing
    end

    return pinv.I
end


function CommonFunctions.write_to_hdf5(pinv::PoincareInvariant2ndCanonical{DT}, sol::Solution, output_file::String) where {DT}
    # h5open(output_file, isfile(output_file) ? "r+" : "w") do h5
    h5open(output_file, "w") do h5
        write(h5, "t", sol.t)
        write(h5, "I", pinv.I)
    end
end


## Common functionality

@generated function compute_integral(Is)
    SU = Ultraspherical(1, 0..1)^2
    Q  = DefiniteIntegral(SU)

    quote
        return Number( $Q * Fun($SU, ApproxFun.transform($SU, Is)) )
    end
end

@generated function compute_derivative1(q)
    Dx = Derivative(Chebyshev(0..1)^2, [1,0])
    CU = Conversion(Ultraspherical(1, 0..1) ⊗ Chebyshev(0..1), Ultraspherical(1, 0..1)^2)

    quote
        return values($CU * ($Dx * q))
    end
end

@generated function compute_derivative2(q)
    Dy = Derivative(Chebyshev(0..1)^2, [0,1])
    CU = Conversion(Chebyshev(0..1) ⊗ Ultraspherical(1, 0..1), Ultraspherical(1, 0..1)^2)

    quote
        return values($CU * ($Dy * q))
    end
end

@generated function compute_canonical_invariant(pinv::Union{PoincareInvariant2nd{DT,ND,NC,NV},PoincareInvariant2ndCanonical{DT,ND,NC,NV}}, q::AbstractArray{DT,2}, p::AbstractArray{DT,2}) where {DT,ND,NC,NV}
    SC = Chebyshev(0..1)^2
    SU = Ultraspherical(1, 0..1)^2

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

@generated function compute_noncanonical_invariant(pinv::PoincareInvariant2nd{DT,ND,NC,NV}, t::DT, q::AbstractArray{DT,2}) where {DT,ND,NC,NV}
    SC = Chebyshev(0..1)^2
    SU = Ultraspherical(1, 0..1)^2
    CU = Conversion(SC, SU)

    Is = zeros(DT, NV)
    qs = zeros(DT, ND, NV)
    qσ = zeros(DT, ND, NV)
    qτ = zeros(DT, ND, NV)
    Ω  = zeros(DT, ND, ND)

    quote
        for j in axes(q,1)
            fq = Fun($SC, ApproxFun.transform($SC, q[j,:]))

            # obtain values of q at the same points as the derivatives
            $qs[j,:] .= itransform($SU, resize!(collect(($CU * fq).coefficients), NC))

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


@generated function compute_noncanonical_correction(pinv::PoincareInvariant2nd{DT,ND,NC,NV}, t::DT, q::AbstractArray{DT,2}, λ::AbstractArray{DT,2}) where {DT,ND,NC,NV}
    SC = Chebyshev(0..1)^2
    SU = Ultraspherical(1, 0..1)^2
    CU = Conversion(SC, SU)

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
            $qs[j,:] .= itransform($SU, resize!(($CU * fq).coefficients, NC))
            $λs[j,:] .= itransform($SU, resize!(($CU * fλ).coefficients, NC))

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
