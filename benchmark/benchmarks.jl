using BenchmarkTools, PoincareInvariants

const SUITE = BenchmarkGroup()

SUITE["FirstPoincareInvariants"] = BenchmarkGroup()

SUITE["SecondPoincareInvariants"] = BenchmarkGroup()

SUITE["SecondPoincareInvariants"]["canonical"] = BenchmarkGroup()
for N in [100 * 100, 1000 * 1000], D in [2, 12, 100]
    name = "SecondPoincareInvariant{$D, Float64}(Ω, $N)"
    pinv = SecondPoincareInvariant{D, Float64}(CanonicalSymplecticMatrix(D), N)
    
    SUITE["SecondPoincareInvariants"]["canonical"][name] =
        @benchmarkable compute($pinv, phasepoints) setup=begin
            phasepoints = rand($(pinv.N), $D)
        end
end

results = run(SUITE)
