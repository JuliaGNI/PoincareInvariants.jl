@safetestset "getpoints, getpointspec and getpointnum" begin
    using PoincareInvariants

    @test getpointspec(1, FirstFinDiffPlan) == 2
    @test getpointspec(50, FirstFinDiffPlan) == 50
    @test getpointspec(345, FirstFinDiffPlan) == 346

    @inferred NTuple{2, Vector{Float32}} getpoints(x -> (x, 2x), Float32, 8, FirstFinDiffPlan)
    @inferred Vector{Float64} getpoints(x -> sin(x), Float64, 124, FirstFinDiffPlan)

    @test getpoints(identity, Float64, 4, FirstFinDiffPlan) ≈ [0, 0.25, 0.5, 0.75]
    @test getpoints(identity, Float64, 6, FirstFinDiffPlan) ≈ [n / 6 for n in 0:5]
    @test getpoints(cos, Float64, 124, FirstFinDiffPlan) ≈ [cos(n / 124) for n in 0:123]

    N = 6342
    testpnts = getpoints(identity, Float64, N, FirstFinDiffPlan)
    pnts3 = getpoints(Float32, N, FirstFinDiffPlan) do x
        x, x^2, sin(x)
    end
    @test pnts3 isa NTuple{3, Vector{Float32}}
    @test all(pnts3) do v
        length(v) == N
    end

    @test pnts3[1] ≈ testpnts
    @test pnts3[2] ≈ testpnts .* testpnts
    @test pnts3[3] ≈ sin.(testpnts)
end

@safetestset "compute! with FirstFinDiffPlan" begin
    using PoincareInvariants
    using StaticArrays: SVector
    using LinearAlgebra: dot

    θ(z, t, p) = SVector{3}(p * z[2], -t * z[1], z[3])
    f(x) = ((s, c) = sincospi(2x); SVector{3}(2c,  5s, c + s))

    N = 6
    fvals = getpoints(f, Float64, N, FirstFinDiffPlan)
    Δx = 1 / N

    testI = sum(1:6) do i
        dfi = ntuple(3) do d
            if i == 1
                return (fvals[d][2] - fvals[d][6]) / (2 * Δx)
            elseif i == 6
                return (fvals[d][1] - fvals[d][5]) / (2 * Δx)
            else
                return (fvals[d][i+1] - fvals[d][i-1]) / (2 * Δx)
            end
        end |> SVector{3}
        w = iseven(i) ? Δx * 4 / 3 : Δx * 2 / 3
        fi = SVector{3}(ntuple(d -> fvals[d][i], 3))
        return w * dot(θ(fi, 2, 3), dfi)
    end

    pinv = FirstPoincareInvariant{Float64, 3}(θ, N, FirstFinDiffPlan())
    @test compute!(pinv, fvals, 2, 3) ≈ testI atol=500eps()
end