@safetestset "getpoints, getpointspec and getpointnum" begin
    using PoincareInvariants

    @test getpointspec(1, FirstFourierPlan) == 1
    @test getpointspec(50, FirstFourierPlan) == 50
    @test getpointspec(345, FirstFourierPlan) == 345

    @inferred Matrix{Float32} getpoints(x -> (x, 2x), Float32, 8, FirstFourierPlan)
    @inferred Vector{Float64} getpoints(x -> sin(x), Float64, 124, FirstFourierPlan)

    @test getpoints(identity, Float64, 4, FirstFourierPlan) ≈ [0, 0.25, 0.5, 0.75]
    @test getpoints(identity, Float64, 6, FirstFourierPlan) ≈ [n / 6 for n in 0:5]
    @test getpoints(cos, Float64, 124, FirstFourierPlan) ≈ [cos(n / 124) for n in 0:123]

    N = 6342
    testpnts = getpoints(identity, Float64, N, FirstFourierPlan)
    pnts3 = getpoints(Float32, N, FirstFourierPlan) do x
        x, x^2, sin(x)
    end
    @test pnts3 isa Matrix{Float32}
    @test size(pnts3) == (N, 3)

    @test pnts3[:, 1] ≈ testpnts
    @test pnts3[:, 2] ≈ testpnts .* testpnts
    @test pnts3[:, 3] ≈ sin.(testpnts)
end

@safetestset "compute! with FirstFourierPlan" begin
    using PoincareInvariants

    T = Float64
    D = 2

    # rotating vector field
    θ(z, p, t) = [-z[2], z[1]]

    # circular path with radius R
    # ḟ = 2π .* R .* (-sinpi(2ϕ), cospi(2ϕ))
    # ∫ ḟ⋅θ dϕ = 2π * R
    R = 0.75
    f(ϕ) = R .* (cospi(2ϕ), sinpi(2ϕ))

    I = 2π * R^2

    N = 10
    plan = FirstFourierPlan{T, D}(θ, N)
    pinv = FirstPoincareInvariant{T, D}(θ, N, plan)
    pnts = getpoints(f, Float64, N, FirstFourierPlan)
    @test I ≈ compute!(pinv, pnts, 0.0, nothing) atol=10eps()
end
