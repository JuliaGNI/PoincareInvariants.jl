
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
