using FastTransforms: PaduaTransformPlan, paduapoints

@testset "PoincareInvariant2 Constructor" begin
    @test PoincareInvariant2 <: AbstractPoincareInvariant

    N = 5
    points = rand(MersenneTwister(321), Float32, 105, 5)
    pinv = PoincareInvariant2(points)

    @test pinv isa PoincareInvariant2{5, Float32, <:PaduaTransformPlan}

    @test pinv.padua_plan isa PaduaTransformPlan
end

@testset "1D Free Particle" begin
    np = 1000
    ppoints = paduapoints(ceil(Int, (-3 + sqrt(1 + 8 * np)) / 2))

    pinv = PoincareInvariant2(ppoints)

    function free_particle!(init, δt)
        init[1] += init[2] .* δt
    end

    Ω(::Any) = [0 -1; 1 0]

    @test compute(pinv, ppoints, Ω) ≈ 4
end

# out = calculate(pinv, phase_points, Ω)
