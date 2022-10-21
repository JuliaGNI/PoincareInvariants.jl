using PoincareInvariants
using CairoMakie  # plotting library

## Pendulum integration ##

pendulum((θ, p)) = (p, -sin(θ))

# struct ForwardEuler end
# timestep((θ, p), dt, ::ForwardEuler) = (θ, p) .+ dt .* pendulum((θ, p))

struct BackwardEuler end
function timestep((θ, p), dt, ::BackwardEuler)
    sinθ, cosθ = sincos(θ)
    p = (p - dt * sinθ) / (1 + dt^2 * cosθ)
    return θ + dt * p, p
end

struct SymplecticEuler end
function timestep((θ, p), dt, ::SymplecticEuler)
    θ = θ + dt * p  # update position
    p = p - dt * sin(θ)  # update momentum
    return θ, p
end

struct RK4 end
function timestep((θ, p), dt, ::RK4)
    k1 = pendulum((θ, p))
    k2 = pendulum((θ, p) .+ 0.5 .* dt .* k1)
    k3 = pendulum((θ, p) .+ 0.5 .* dt .* k2)
    k4 = pendulum((θ, p) .+ dt .* k3)
    return (θ, p) .+ (dt/6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

"""
    integrate((θ₀, p₀), dt, nsteps, nt, method)

start at `(θ₀, p₀)` and integrate the equations of motion using `method`.
Returns the timeseries as a vector of tuples. `nt` points are saved,
`nsteps` steps are taken from saved point to saved point and `dt` is the
size of each time step.
"""
function integrate((θ₀, p₀), dt, nsteps, nt, method)
    out = Vector{Tuple{Float64, Float64}}(undef, nt)
    (θ, p) = out[1] = (θ₀, p₀)
    for i in 2:nt
        for _ in 1:nsteps
            (θ, p) = timestep((θ, p), dt, method)
        end
        out[i] = (θ, p)
    end
    return out
end

integrate(mat::AbstractMatrix, dt, nsteps, nt, method) = map(eachrow(mat)) do pnt
    integrate((pnt[1], pnt[2]), dt, nsteps, nt, method)
end

## Invariants ##

pi1 = CanonicalFirstPI{Float64, 2}(500)
pi2 = CanonicalSecondPI{Float64, 2}(10_000)

@assert 500 ≤ getpointnum(pi1) ≤ 750
@assert 10_000 ≤ getpointnum(pi2) ≤ 15_000

I1 = 3π
pnts1 = getpoints(pi1) do ϕ
    sinpi(2ϕ), 3 * cospi(2ϕ)
end

I2 = 16
pnts2 = getpoints(pi2) do x, y
    4 .* (x, y) .- 2
end

# get points used in domain
getpoints(pi1)
getpoints(pi2)

# plot initial points in phase space
let fig = Figure(resolution=(450, 700))
    ax1 = Axis(fig[1, 1], xlabel="θ", ylabel="p", aspect=DataAspect())
    scatter!(ax1, pnts1[:, 1], pnts1[:, 2]; markersize=2, label="curve")
    scatter!(ax1, pnts2[:, 1], pnts2[:, 2]; markersize=2, label="surface")
    axislegend(ax1; position=:lt)
    save("pendulum_init.png", fig)
    fig
end

@assert isapprox(compute!(pi1, pnts1), I1; atol=10eps())
@assert isapprox(compute!(pi2, pnts2), I2; atol=10eps())

series1 = integrate(pnts1, 0.05, 15, 5, SymplecticEuler())
@assert all(compute!(pi1, series1)) do I
    abs(I - I1) < 10^(-14)
end

series2 = integrate(pnts2, 0.05, 15, 5, SymplecticEuler())
@assert all(compute!(pi2, series2)) do I
    abs(I - I2) < 10^(-13)
end

calcerrs(Is, I0) = max.(abs.(Is .- I0) ./ eps(), 1.0) .* eps()

let fig = Figure(resolution=(2*800, 2*400)), dt = 0.05, nsteps = 15, nt = 5
    times = range(0.0; step=dt * nsteps, length=nt)

    azimuth = π/2 - 0.05
    xlabel = "θ"
    ylabel = "p"
    zlabel = "t"
    ax1 = Axis3(fig[1, 1]; elevation=0.22, azimuth=azimuth, aspect=(2, 1, 1),
        xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)
    ax2 = Axis3(fig[1, 2]; elevation=0.22, azimuth=azimuth, aspect=(2, 1, 1),
        xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)

    Iax1 = Axis(fig[2, 1]; yscale=log10, xlabel="time", ylabel="Absolute Error (log scale)")
    Iax2 = Axis(fig[2, 2]; yscale=log10, xlabel="time", ylabel="Absolute Error (log scale)")

    for (i, method) in enumerate([SymplecticEuler(), RK4(), BackwardEuler()])
        intdata1 = integrate(getpoints(ϕ -> (sinpi(2ϕ), 3 * cospi(2ϕ)), pi1), dt, nsteps, nt, method)
        intdata2 = integrate(getpoints((x, y) -> 4 .* (x, y) .- 2, pi2), dt, nsteps, nt, method)

        for (j, t) in enumerate(times)
            pnts1x = map(x -> x[j][1], intdata1)
            pnts1y = map(x -> x[j][2], intdata1)
            zs1 = fill(t, getpointnum(pi1))
            lines!(ax1, pnts1x, pnts1y, zs1; color=Cycled(i))

            zs2 = fill(t, getpointnum(pi2))
            pnts2x = map(x -> x[j][1], intdata2)
            pnts2y = map(x -> x[j][2], intdata2)
            scatter!(ax2, pnts2x, pnts2y, zs2; color=Cycled(i), markersize=0.5)
        end

        lines!(Iax1, times, calcerrs(compute!(pi1, intdata1), I1); color=Cycled(i), label=string(method)[1:end-2])
        lines!(Iax2, times, calcerrs(compute!(pi2, intdata2), I2); color=Cycled(i), label=string(method)[1:end-2])
    end

    axislegend(Iax2; position=:lt)

    save("pendulum.png", fig)
    fig
end
