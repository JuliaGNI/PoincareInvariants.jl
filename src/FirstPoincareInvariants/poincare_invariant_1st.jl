
struct PoincareInvariant1st{ET,DT,TT,ΘT}
    equ::ET
    Θ::ΘT
    Δt::TT
    nloop::Int
    ntime::Int
    nsave::Int
    nt::Int
    I::OffsetArray{DT,1,Vector{DT}}
    J::OffsetArray{DT,1,Vector{DT}}
    L::OffsetArray{DT,1,Vector{DT}}
end

function PoincareInvariant1st(f_equ::Function, f_loop::Function, Θ::ΘT, Δt::TT, d::Int, nloop::Int, ntime::Int, nsave::Int=1, DT=Float64) where {TT,ΘT}

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
    q₀ = [f_loop(i/nloop) for i in 1:nloop]

    # initialise euation
    equ = f_equ(q₀)

    # create arrays for results
    nt = div(ntime, nsave)

    I = OffsetArray(zeros(DT, nt+1), 0:nt)
    J = OffsetArray(zeros(DT, nt+1), 0:nt)
    L = OffsetArray(zeros(DT, nt+1), 0:nt)

    PoincareInvariant1st{typeof(equ),DT,TT,ΘT}(equ, Θ, Δt, nloop, ntime, nsave, nt, I, J, L)
end


function evaluate_poincare_invariant(pinv::PoincareInvariant1st, sol::Solution)
    p = zeros(size(sol.q[begin],1), size(sol.q,1), size(sol.q,2))
    g = zeros(size(sol.q[begin],1), size(sol.q,1), size(sol.q,2))
    v = zeros(size(sol.q[begin],1), size(sol.q,2))
    γ = zeros(size(sol.q[begin],1), size(sol.q,2))

    compute_one_form(sol.t, sol.q, p, pinv.Θ)

    if isdefined(sol, :λ)
        compute_correction(sol.t, sol.q, sol.λ, g, pinv.equ.g)
    end

    for i in axes(sol.q,1)
        compute_velocity(hcat(sol.q[i,:]...), v)
        pinv.I[i] = compute_loop_integral(hcat(sol.p[i,:]...), v)

        if isdefined(sol, :p)
            pinv.J[i] = compute_loop_integral(hcat(sol.d[i,:]...), v)
        end

        if isdefined(sol, :λ)
            compute_velocity(hcat(sol.λ[i,:]...), γ)
            pinv.L[i] = compute_loop_integral(p[:,i,:] .- pinv.Δt .* g[:,i,:], v .- pinv.Δt .* γ)
        end
    end

    (pinv.I, pinv.J, pinv.L)
end


function write_to_hdf5(pinv::PoincareInvariant1st, sol::Solution, output_file::String)
    # h5open(output_file, isfile(output_file) ? "r+" : "w") do h5
    h5open(output_file, "w") do h5

        write(h5, "t", sol.t)
        write(h5, "I", pinv.I)

        isdefined(sol, :p) ? write(h5, "J", pinv.J) : nothing
        isdefined(sol, :λ) ? write(h5, "L", pinv.L) : nothing

    end
end
