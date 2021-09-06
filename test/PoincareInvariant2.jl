using FastTransforms: PaduaTransformPlan

@testset "PoincareInvariant2" begin
    @test PoincareInvariant2 <: AbstractPoincareInvariant

    padua_points = get_padua_points(100)

    f(v) = SVector{2, Float32}(cos(v[1]) * sin(v[2]), v[1] + v[2])

    phase_points = f.(padua_points)

    tup = (
        first.(phase_points),
        last.(phase_points)
    )

    pinv = PoincareInvariant2(tup)

    @test pinv isa PoincareInvariant2{2, Float32, <:PaduaTransformPlan}
end

# out = calculate(pinv, phase_points, Ω)
