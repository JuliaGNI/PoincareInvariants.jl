
module PoincareInvariant2ndTest

    using GeometricIntegrators
    using PoincareInvariants
    using SymPy

    const Δt = 10.
    const B₀ = 1.
    const r₀ = 0.5
    const z₀ = 0.0
    const z₁ = 0.1
    const u₀ = 5E-1
    const u₁ = 5E-2
    const kz = 2
    const ku = 4


    function B(t, q)
        B₀ * (1 + q[1]^2 + q[2]^2)
    end

    function gcϑ(t, q, ϑ)
        ϑ[1] =-q[2] * (1 + B(t,q)) / 4
        ϑ[2] =+q[1] * (1 + B(t,q)) / 4
        ϑ[3] = q[4]
        ϑ[4] = 0
        nothing
    end

    function gcϑ(t, q)
        ϑ = zero(q)
        gcϑ(t, q, ϑ)
        return ϑ
    end

    function gcω(t, q, Ω)
        Ω[1,1] = 0
        Ω[1,2] =-B(t,q)
        Ω[1,3] = 0
        Ω[1,4] = 0

        Ω[2,1] =+B(t,q)
        Ω[2,2] = 0
        Ω[2,3] = 0
        Ω[2,4] = 0

        Ω[3,1] = 0
        Ω[3,2] = 0
        Ω[3,3] = 0
        Ω[3,4] =+q[4]

        Ω[4,1] = 0
        Ω[4,2] = 0
        Ω[4,3] =-q[4]
        Ω[4,4] = 0

        nothing
    end

    function gc_surface_q(s,t)
        x  = r₀*(s-0.5)
        y  = r₀*(t-0.5)
        z  = z₀ + z₁ * cos(2π*kz*s) * cos(2π*kz*t)
        u  = u₀ + u₁ * sin(2π*ku*s) * sin(2π*ku*t)

        [x, y, z, u]
    end

    function gc_surface_p(s,t)
        gcϑ(0, gc_surface_q(s,t))
    end

    function gc_dummy(t, a, b, c)
        nothing
    end

    function gc_dummy_pode(q₀, p₀)
        PODE(gc_dummy, gc_dummy, q₀, p₀)
    end

    function gc_dummy_iode(q₀)
        p₀ = [gcϑ(0,q) for q in q₀]
        IODE(gc_dummy, gc_dummy, gc_dummy, q₀, p₀; v̄=gc_dummy)
    end


    function compute_canonical_invariant_approxfun(nx, ny, nt=0)
        pinv = PoincareInvariant2ndCanonical(gc_dummy_pode, gc_surface_q, gc_surface_p, Δt, 4, nx, ny, nt)
        sol  = Solution(pinv.equ, Δt, nt)

        I = evaluate_poincare_invariant(pinv, sol)

        return I[0]
    end


    function compute_invariant_approxfun(nx, ny, nt=0)
        pinv = PoincareInvariant2ndApproxFun(gc_dummy_iode, gc_surface_q, gcω, Δt, 4, nx, ny, nt)
        sol  = Solution(pinv.equ, Δt, nt)

        I, J, ΔI, ΔJ = evaluate_poincare_invariant(pinv, sol)

        return I[0], J[0]
    end


    function compute_invariant_opq(nx, ny, nt=0)
        pinv = PoincareInvariant2ndOPQ(gc_dummy_iode, gc_surface_q, gcω, Δt, 4, nx, ny, nt)
        sol  = Solution(pinv.equ, Δt, nt)

        I, J, ΔI, ΔJ = evaluate_poincare_invariant(pinv, sol)

        return I[0], J[0]
    end


    function compute_invariant_trapezoidal(nx, ny, nt=0)
        pinv = PoincareInvariant2ndTrapezoidal(gc_dummy_iode, gc_surface_q, gcω, Δt, 4, nx, ny, nt)
        sol  = Solution(pinv.equ, Δt, nt)

        I, J, ΔI, ΔJ = evaluate_poincare_invariant(pinv, sol)

        return I[0], J[0]
    end


    function compute_invariant_analytical()
        local verbosity = get_config(:verbosity)

        if get_config(:verbosity) > 1
            println()
            println("Second Euler-Poincaré Integral Invariant (Analytical)")
            println("=====================================================")
            println()
        end
        
        s, t = Sym("s, t")
        x, y, z, u = gc_surface_q(s,t)

        S = simplify(B₀ * (1 + x^2 + y^2) * ( diff(x, s) * diff(y, t) + diff(x, t) * diff(y, s) )
                   + u * ( diff(z, s) * diff(u, t) + diff(z, t) * diff(u, s) ) )
        I = N(SymPy.integrate(S, (s, 0, 1), (t, 0, 1)), 16)

        verbosity > 1 ? println("           I = ", I) : nothing

        return I
    end


    export compute_canonical_invariant_approxfun,
           compute_invariant_approxfun,
           compute_invariant_opq,
           compute_invariant_trapezoidal,
           compute_invariant_analytical

end
