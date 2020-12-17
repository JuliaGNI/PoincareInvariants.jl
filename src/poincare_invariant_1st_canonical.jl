
struct PoincareInvariant1stCanonical{ET,DT,TT}
    equ::ET
    Δt::TT
    nloop::Int
    ntime::Int
    nsave::Int
    nt::Int
    I::OffsetArray{DT,1,Vector{DT}}
end

function PoincareInvariant1stCanonical(f_equ::Function, f_loop_q::Function, f_loop_p::Function, Δt::TT, d::Int, nloop::Int, ntime::Int, nsave::Int, DT=Float64) where {TT}

    if get_config(:verbosity) > 1
        println()
        println("First Euler-Poincaré Integral Invariant")
        println("=======================================")
        println()
        println(" nloop = ", nloop)
        println(" ntime = ", ntime)
        println(" nsave = ", nsave)
        println(" Δt    = ", Δt)
        println()
    end

    # compute initial conditions
    q₀ = [f_loop_q(i/nloop) for i in 1:nloop]
    p₀ = [f_loop_p(i/nloop) for i in 1:nloop]

    # initialise euation
    equ = f_equ(q₀, p₀)

    # create arrays for results
    nt = div(ntime, nsave)

    I = OffsetArray(zeros(DT, nt+1), 0:nt)

    PoincareInvariant1stCanonical{typeof(equ),DT,TT}(equ, Δt, nloop, ntime, nsave, nt, I)
end


function evaluate_poincare_invariant(pinv::PoincareInvariant1stCanonical, sol::SolutionPODE)
    v = zeros(size(sol.q[begin],1), size(sol.q,2))

    for i in axes(sol.q,1)
        compute_velocity(hcat(sol.q[i,:]...), v)
        pinv.I[i] = compute_loop_integral(hcat(sol.p[i,:]...), v)
    end

    pinv.I
end


function write_to_hdf5(pinv::PoincareInvariant1stCanonical, sol::Solution, output_file::String)
    # h5open(output_file, isfile(output_file) ? "r+" : "w") do h5
    h5open(output_file, "w") do h5
        write(h5, "t", sol.t)
        write(h5, "I", pinv.I)
    end
end
