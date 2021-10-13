using BenchmarkTools, PoincareInvariants

const SUITE = BenchmarkGroup()

SUITE["FirstPoincareInvariants"] = BenchmarkGroup()

SUITE["SecondPoincareInvariants"] = BenchmarkGroup()

SUITE["SecondPoincareInvariants"]["canonical"] = BenchmarkGroup()
for N in [100 * 100, 1000 * 1000], D in [2, 12, 100]
    name = "SecondPoincareInvariant{Float64}(Ω, $D, $N)"

    SUITE["SecondPoincareInvariants"]["canonical"][name] =
        @benchmarkable compute!(pinv, phasepoints, 0, nothing) setup=begin
            Ω(z, t, p) = CanonicalSymplecticMatrix($D)
            pinv = SecondPoincareInvariant{Float64}(Ω, $D, $N, Val(false))
            phasepoints = rand(getpointnum(pinv), $D)
        end
end

results = run(SUITE)
