using PoincareInvariants
using StaticArrays
using CairoMakie

E(x, y) = (-x, -y^3)

A(x, y) = 5 .* (-y, x)
B(x, y) = 10.0
Bx(x, y, Δx) = 10 * Δx  # integral in x direction from x to x + Δx
By(x, y, Δy) = 10 * Δy  # integral in y direction from y to y + Δy

struct Split2 end

function ϕx((x, y, vx, vy), dt)
    Δx = vx * dt
    return (x + Δx, y, vx, vy - Bx(x, y, Δx))
end

function ϕy((x, y, vx, vy), dt)
    Δy = vy * dt
    return (x, y + Δy, vx + By(x, y, Δy), vy)
end

function ϕE((x, y, vx, vy), dt)
    (Δvx, Δvy) = dt .* E(x, y)
    return (x, y, vx + Δvx, vy + Δvy)
end

function timestep(z, dt, ::Split2)
    hdt = 0.5 * dt
    z = ϕx(z, hdt)
    z = ϕy(z, hdt)
    z = ϕE(z, dt)
    z = ϕy(z, hdt)
    z = ϕx(z, hdt)
    return z
end

function zdot((x, y, vx, vy))
    ex, ey = E(x, y); b = B(x, y)
    (vx, vy, ex + b * vy, ey - b * vx)
end

struct Euler end
timestep(z, dt, ::Euler) = z .+ dt .* zdot(z)

struct RK4 end
function timestep(z, dt, ::RK4)
    hdt = 0.5 .* dt
    k1 = zdot(z)
    k2 = zdot(z .+ hdt .* k1)
    k3 = zdot(z .+ hdt .* k2)
    k4 = zdot(z .+  dt .* k3)
    return z .+ dt .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4) .* (1/6)
end

"""
    integrate(z0, dt, nsteps, nt, method)

start at `z0` and integrate the equations of motion using `method`.
Returns the timeseries as a vector of tuples. `nt` points are saved,
`nsteps` steps are taken from saved point to saved point and `dt` is the
size of each time step.
"""
function integrate(z0, dt, nsteps, nt, method)
    out = Vector{NTuple{4, Float64}}(undef, nt)
    z = out[1] = z0
    for i in 2:nt
        for _ in 1:nsteps
            z = timestep(z, dt, method)
        end
        out[i] = z
    end
    return out
end

integrate(mat::AbstractMatrix, dt, nsteps, nt, method) = map(eachrow(mat)) do r
    integrate((r[1], r[2], r[3], r[4]), dt, nsteps, nt, method)
end

## Invariants ##

function oneform((x, y, vx, vy), ::Real, ::Any)
    p = (vx, vy) .+ A(x, y)
    @SVector [p[1], p[2], 0, 0]
end

function twoform(z, ::Real, ::Any)
    b = B(z[1], z[2])
    @SMatrix [ 0  b -1  0;
              -b  0  0 -1;
               1  0  0  0;
               0  1  0  0]
end

pi1 = FirstPI{Float64, 4}(oneform, 1_000)
pi2 = SecondPI{Float64, 4}(twoform, 10_000)

I1 = 0
pnts1 = getpoints(pi1) do θ
    10 .* (0, 0, sinpi(2θ), cospi(2θ))
end

I2 = 1_000
pnts2 = getpoints(pi2) do x, y
    10 .* (y - 0.5, x - 0.5, 0, 0)
end

@assert isapprox(compute!(pi1, pnts1, 52.3, "optional parameter"), I1; atol=10^(-15))
@assert isapprox(compute!(pi2, pnts2), I2; atol=10^(-11))

times = range(0.0; step=0.05 * 50, length=5)

series1 = integrate(pnts1, 0.05, 50, 5, Split2())
@assert all(compute!(pi1, series1, times, ("optional parameters", 3.7, 42))) do I
    abs(I - I1) < 10^(-13)
end

series2 = integrate(pnts2, 0.05, 50, 5, Split2())
@assert all(compute!(pi2, series2)) do I
    abs(I - I2) < 10^(-11)
end

calcerrs(Is, I0) = max.(abs.(Is .- I0) ./ eps(), 1.0) .* eps()

let fig = Figure(resolution=(2*800, 2*400)), dt = 0.05, nsteps = 50, nt = 5
    times = range(0.0; step=dt * nsteps, length=nt)

    azimuth = π/2 - 0.05
    xlabel = "x"
    ylabel = "y"
    zlabel = "t"
    ax1 = Axis3(fig[1, 1]; elevation=0.22, azimuth=azimuth, aspect=(2, 1, 1),
        xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)
    ax2 = Axis3(fig[1, 2]; elevation=0.22, azimuth=azimuth, aspect=(2, 1, 1),
        xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)

    Iax1 = Axis(fig[2, 1]; yscale=log10, xlabel="time", ylabel="Absolute Error (log scale)")
    Iax2 = Axis(fig[2, 2]; yscale=log10, xlabel="time", ylabel="Absolute Error (log scale)")

    for (i, method) in enumerate([RK4(), Split2()])
        intdata1 = integrate(getpoints(θ -> 10 .* (0, 0, sinpi(2θ), cospi(2θ)), pi1), dt, nsteps, nt, method)
        intdata2 = integrate(getpoints((x, y) -> 10 .* (y - 0.5, x - 0.5, 0, 0), pi2), dt, nsteps, nt, method)

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

    save("particle.png", fig)
    fig
end
