using PoincareInvariants
using CairoMakie
using StatsBase: Histogram

## Maps ##

function wavy(x, y)
    x, y = x, y + 2 * sinpi(x)  # make wavy
    x, y = x - y * exp(-y^2 / 4), y
    x, y = -y, x  # rotate
    return x, y
end

function wavy(init, pnt, ::Val{N}) where N
    x, y = init(pnt...)
    for _ in 1:N
        x, y = wavy(x, y)
    end
    return x, y
end

calcerrs(Is, I0) = max.(abs.(Is .- I0) ./ eps(), 1.0) .* eps()

# Makie doesn't handle plotting a million+ points very well, so I make them into a heatmap
function makehist(f, N, edges)
    h = Histogram(edges)
    n = ceil(Int, sqrt(N))
    xs = range(0, 1, length=n)
    for x in xs
        for y in xs
            push!(h, f(x, y))
        end
    end
    return h
end

init2(x, y) = 2 .* (x, y) .- 1

let N = 3, fig = Figure(resolution=(4*800, (N+1)*800)), I0 = 4
    pinvnums = [100, 1_000, 10_000, 100_000, 1_000_000]
    pinvs = CanonicalSecondPI{Float64, 2}.(pinvnums)

    function calcI2(n, pinv)
        f(x, y) = wavy(init2, (x, y), Val(n))
        compute!(pinv, getpoints(f, pinv))
    end

    Is = [calcI2(n, pinv) for n in 0:N, pinv in pinvs]
    errs = calcerrs(Is, I0)

    xnums = [100, 100, 100, 5000, 1000]
    msizes = [7.5, 5, 2, 1.5]
    for n in 0:N
        f(x, y) = wavy(init2, (x, y), Val(n))
        h = makehist(f, 1000_000, (-5:0.01:5, -5:0.01:5))

        hax = Axis(fig[n+2, 2])
        plot!(hax, h; colormap=:grayC, lowclip=RGBAf(0,0,0,0), colorrange=(eps(), 20))

        xs = range(0, 1, length=xnums[n+1])
        zs = [f(x, y) for x in xs, y in xs]
        fxs = first.(zs)
        fys = last.(zs)

        xax = Axis(fig[n+2, 3])
        yax = Axis(fig[n+2, 4])
        heatmap!(xax, xs, xs, fxs; colormap=:viridis, colorrange=(-5, 5))
        heatmap!(yax, xs, xs, fys; colormap=:viridis, colorrange=(-5, 5))
        if n == 0
            ms = 10
            scatter!(xax, getpoints(pinvs[n+1]); color=:white, markersize=ms)
            scatter!(yax, getpoints(pinvs[n+1]); color=:white, markersize=ms)
        end

        errax = Axis(fig[n+2, 5]; xscale=log10, yscale=log10,
            xlabel="Number of points (log scale)",
            ylabel="Absolute Error (log scale)")
        n ≥ 1 && linkyaxes!(errax, content(fig[n+1, 5]))
        scatter!(errax, pinvnums, errs[n+1, :]; color=:black, markersize=10)
    end

    txsz = 26
    Label(fig[1, 2], "surface"; textsize = txsz, tellwidth=false)
    Label(fig[1, 3], "x component of parameterisation"; textsize = txsz, tellwidth=false)
    Label(fig[1, 4], "y component of parameterisation"; textsize = txsz, tellwidth=false)
    Label(fig[1, 5], "surface area error"; textsize = txsz, tellwidth=false)

    Label(fig[2, 1], "undeformed"; rotation = pi/2, textsize = txsz, tellheight=false)
    Label(fig[3, 1], "once wavy"; rotation = pi/2, textsize = txsz, tellheight=false)
    Label(fig[4, 1], "twice wavy"; rotation = pi/2, textsize = txsz, tellheight=false)
    Label(fig[5, 1], "thrice wavy"; rotation = pi/2, textsize = txsz, tellheight=false)

    save("wavy2.png", fig)
    nothing
end

init1(θ) = (cospi(2θ), -sinpi(2θ))

let N = 4, fig = Figure(resolution=(4*800, (N+1)*800)), I0 = π
    pinvnums = [100, 1_000, 10_000, 100_000]
    pinvs = CanonicalFirstPI{Float64, 2}.(pinvnums)

    function calcI1(n, pinv)
        f(θ) = wavy(init1, θ, Val(n))
        compute!(pinv, getpoints(f, pinv))
    end

    Is = [calcI1(n, pinv) for n in 0:N, pinv in pinvs]
    errs = calcerrs(Is, I0)

    for n in 0:N
        f(θ) = wavy(init1, θ, Val(n))
        θs = getpoints(pinvs[end])
        pnts = getpoints(f, pinvs[end])

        linax = Axis(fig[n+2, 2]; limits=(-5, 5, -5, 5))
        lines!(linax, pnts[:, 1], pnts[:, 2])

        pax = Axis(fig[n+2, 3:4])
        lines!(pax, θs, pnts[:, 1]; label="x component")
        lines!(pax, θs, pnts[:, 2]; label="y component")
        axislegend(pax; position=:rb)

        errax = Axis(fig[n+2, 5]; xscale=log10, yscale=log10,
            xlabel="Number of points (log scale)",
            ylabel="Absolute Error (log scale)")
        n ≥ 1 && linkyaxes!(errax, content(fig[n+1, 5]))
        scatter!(errax, pinvnums, errs[n+1, :]; color=:black, markersize=10)
    end

    txsz = 26
    Label(fig[1, 2], "loop"; textsize = txsz, tellwidth=false)
    Label(fig[1, 3:4], "parameterisation"; textsize = txsz, tellwidth=false)
    Label(fig[1, 5], "first invariant error"; textsize = txsz, tellwidth=false)

    Label(fig[2, 1], "undeformed"; rotation = pi/2, textsize = txsz, tellheight=false)
    Label(fig[3, 1], "once wavy"; rotation = pi/2, textsize = txsz, tellheight=false)
    Label(fig[4, 1], "twice wavy"; rotation = pi/2, textsize = txsz, tellheight=false)
    Label(fig[5, 1], "thrice wavy"; rotation = pi/2, textsize = txsz, tellheight=false)
    Label(fig[6, 1], "four times wavy"; rotation = pi/2, textsize = txsz, tellheight=false)

    save("wavy1.png", fig)
    nothing
end
