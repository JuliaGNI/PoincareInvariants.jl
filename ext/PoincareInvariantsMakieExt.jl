module PoincareInvariantsMakieExt

using Makie
using PoincareInvariants
using PoincareInvariants: FirstPoincareInvariant, SecondPoincareInvariant,
    AbstractPoincareInvariant, compute!, getpointspec
using GeometricSolutions: EnsembleSolution, nsamples, ntime

# shared 3D view for the advected curve/surface plots
const VIEW3 = (; azimuth = 1.775π, elevation = π / 8, aspect = (1, 1, 1.6))

saved_times(sol) = [sol[1].t[n] for n in 0:ntime(sol[1])]
relerr(Is) = (Is .- Is[1]) ./ Is[1]

# `nsteps` evenly-spaced saved time indices (0-based)
timeslices(sol, nsteps) = round.(Int, range(0, ntime(sol[1]), length = nsteps))

# 2D phase space points of every ensemble member at saved time index n; works for both
# canonical (q, p) solutions and noncanonical solutions where the phase point is q
function phasepoints(sol, n)
    m1 = sol[1]
    if hasproperty(m1, :p) && length(m1.q[n]) == 1
        ([sol[j].q[n][1] for j in 1:nsamples(sol)], [sol[j].p[n][1] for j in 1:nsamples(sol)])
    else
        ([sol[j].q[n][1] for j in 1:nsamples(sol)], [sol[j].q[n][2] for j in 1:nsamples(sol)])
    end
end

# symbol used in the axis label of the relative-error plot
invariant_symbol(::FirstPoincareInvariant) = "I₁"
invariant_symbol(::SecondPoincareInvariant) = "I₂"

## invariant error traces ##########################################################

function PoincareInvariants.plot_invariant!(ax, pinv::AbstractPoincareInvariant,
        sol::EnsembleSolution; p = nothing, label = nothing)
    ts = saved_times(sol)
    Is = compute!(pinv, sol, p)
    scatter!(ax, ts, relerr(Is); label = label)
    Makie.xlims!(ax, first(ts), last(ts))
    ax
end

function _invariant_axis(fig, pinv, title)
    s = invariant_symbol(pinv)
    ax = Axis(fig[1, 1]; xlabel = "t",
        ylabel = "Relative Error ($(s)(t)-$(s)(0))/$(s)(0)", title = title)
    hlines!(ax, [0.0]; color = :gray, linestyle = :dash)
    ax
end

function PoincareInvariants.plot_invariant(pinv::AbstractPoincareInvariant,
        sol::EnsembleSolution; p = nothing, title = "")
    fig = Figure()
    ax = _invariant_axis(fig, pinv, title)
    plot_invariant!(ax, pinv, sol; p = p)
    fig
end

function PoincareInvariants.plot_invariant(pinv::AbstractPoincareInvariant,
        sols::Pair...; p = nothing, title = "", position = :lb)
    fig = Figure()
    ax = _invariant_axis(fig, pinv, title)
    for (label, sol) in sols
        plot_invariant!(ax, pinv, sol; p = p, label = label)
    end
    axislegend(ax; position = position)
    fig
end

## advected loop ###################################################################

function PoincareInvariants.plot_loop!(ax, sol::EnsembleSolution; nsteps = 10)
    ts = saved_times(sol)
    pal = Makie.wong_colors()
    for (k, n) in enumerate(timeslices(sol, nsteps))
        xs, ys = phasepoints(sol, n)
        lines!(ax, xs, ys, fill(ts[n + 1], length(xs)); color = pal[mod1(k, length(pal))])
    end
    Makie.zlims!(ax, first(ts), last(ts))
    ax
end

function PoincareInvariants.plot_loop(sol::EnsembleSolution;
        xlabel = "q₁", ylabel = "q₂", title = "Advected Loop", nsteps = 10)
    fig = Figure()
    ax = Axis3(fig[1, 1]; xlabel = xlabel, ylabel = ylabel, zlabel = "t", VIEW3..., title = title)
    plot_loop!(ax, sol; nsteps = nsteps)
    fig
end

## advected surface ################################################################

function PoincareInvariants.plot_surface!(ax, grid::SecondPoincareInvariant,
        sol::EnsembleSolution; nsteps = 10)
    nx, ny = getpointspec(grid)
    ts = saved_times(sol)
    pal = Makie.wong_colors()
    for (k, n) in enumerate(timeslices(sol, nsteps))
        xs, ys = phasepoints(sol, n)
        surface!(ax, reshape(xs, ny, nx), reshape(ys, ny, nx), fill(ts[n + 1], ny, nx);
            color = fill(pal[mod1(k, length(pal))], ny, nx), shading = NoShading)
    end
    Makie.zlims!(ax, first(ts), last(ts))
    ax
end

function PoincareInvariants.plot_surface(grid::SecondPoincareInvariant, sol::EnsembleSolution;
        xlabel = "q₁", ylabel = "q₂", title = "Advected Surface", nsteps = 10)
    fig = Figure()
    ax = Axis3(fig[1, 1]; xlabel = xlabel, ylabel = ylabel, zlabel = "t", VIEW3..., title = title)
    plot_surface!(ax, grid, sol; nsteps = nsteps)
    fig
end

end  # module
