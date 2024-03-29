using BenchmarkTools, PoincareInvariants

const SUITE = BenchmarkGroup()

SUITE["FirstPoincareInvariants"] = BenchmarkGroup()

SUITE["SecondPoincareInvariants"] = BenchmarkGroup()

SUITE["SecondPoincareInvariants"]["canonical"] = BenchmarkGroup()
for N in [100 * 100, 1000 * 1000], D in [2, 12, 100]
    name = "SecondPoincareInvariant{Float64}(ω, $D, $N)"

    SUITE["SecondPoincareInvariants"]["canonical"][name] =
        @benchmarkable compute!(pinv, phasepoints, 0, nothing) setup=begin
            ω(z, t, p) = canonical_two_form($D)
            pinv = SecondPoincareInvariant{Float64}(ω, $D, $N)
            phasepoints = getpoints(pinv) do x, y
                rand($D)
            end
        end
end

results = run(SUITE)
