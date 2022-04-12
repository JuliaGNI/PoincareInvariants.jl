using PoincareInvariants.SecondFinDiffPlans

@safetestset "getpoints, getpointspec and getpointnum" begin
    using PoincareInvariants

    @test getpointspec(50, SecondFinDiffPlan) == (9, 9)
    @test getpointspec(87, SecondFinDiffPlan) == (11, 11)
    @test getpointspec((53, 42), SecondFinDiffPlan) == (53, 43)

    @test getpointnum((3, 5)) == 15
    @test getpointnum((423789, 326)) == 423789 * 326
    @test getpointnum((123, 908343)) == 123 * 908343

    @inferred Matrix{Float64} getpoints((x, y) -> (x, x+y), Float64, (7, 5), SecondFinDiffPlan)
    @inferred Vector{Float32} getpoints((x, y) -> sin(x), Float32, (9, 9), SecondFinDiffPlan)

    @test getpoints((x, y) -> (x, y), Float64, (5, 3), SecondFinDiffPlan) ≈ [
        0    0  ;
        0    0.5;
        0    1  ;
        0.25 0  ;
        0.25 0.5;
        0.25 1  ;
        0.5  0  ;
        0.5  0.5;
        0.5  1  ;
        0.75 0  ;
        0.75 0.5;
        0.75 1  ;
        1.0  0  ;
        1.0  0.5;
        1.0  1  ]

    @test getpoints((x, y) -> (4x, 6y), Float64, (5, 7), SecondFinDiffPlan) ≈ [
        0 0;
        0 1;
        0 2;
        0 3;
        0 4;
        0 5;
        0 6;
        1 0;
        1 1;
        1 2;
        1 3;
        1 4;
        1 5;
        1 6;
        2 0;
        2 1;
        2 2;
        2 3;
        2 4;
        2 5;
        2 6;
        3 0;
        3 1;
        3 2;
        3 3;
        3 4;
        3 5;
        3 6;
        4 0;
        4 1;
        4 2;
        4 3;
        4 4;
        4 5;
        4 6]

    @test getpoints((x, y) -> y, Float64, (11, 11), SecondFinDiffPlan) ≈ Float64[
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    ]

    ps = getpointspec(6342, SecondFinDiffPlan)
    testpnts = getpoints((x, y) -> (x, y), Float64, ps, SecondFinDiffPlan)
    pnts4 = getpoints(Float64, ps, SecondFinDiffPlan) do x, y
        x, y, x * y, x + y
    end
    @test pnts4 isa Matrix{Float64}
    @test size(pnts4) == (getpointnum(ps), 4)

    @test pnts4[:, 1] ≈ testpnts[:, 1]
    @test pnts4[:, 2] ≈ testpnts[:, 2]
    @test pnts4[:, 3] ≈ testpnts[:, 1] .* testpnts[:, 2]
    @test pnts4[:, 4] ≈ testpnts[:, 1] .+ testpnts[:, 2]
end

@safetestset "Differentiation" begin
    using ..SecondFinDiffPlans: getpoints, SecondFinDiffPlan, differentiate

    # Finite difference derivatives should be exact for a quadratic
    let f(x, y)  = 1 + 0.75*x + y - 0.5*x*x - y*y - 1.5*x*y
        fx(x, y) = 0.75 - x - 1.5*y
        fy(x, y) = 1 - 2*y - 1.5*x

        for (nx, ny, maxerr) in [(3, 3, 50), (21, 33, 500), (331, 123, 5000)]
            vals = getpoints(f, Float64, (nx, ny), SecondFinDiffPlan)

            @test all((ix, iy) for ix in 1:nx, iy in 1:ny) do (ix, iy)
                ∂x, ∂y = differentiate(vals, ix, iy, (nx, ny))
                x ,  y = (ix - 1) / (nx - 1), (iy - 1) / (ny - 1)
                tx, ty = fx(x, y), fy(x, y)
                return abs(tx - ∂x) / eps() < maxerr && abs(ty - ∂y) / eps() < maxerr
            end
        end
    end

    let f(x, y) = sinpi(x) * exp(y)
        fx(x, y) = π * cospi(x) * exp(y)
        fy(x, y) = sinpi(x) * exp(y)

        nx, ny = 1321, 1235
        vals = getpoints(f, Float64, (nx, ny), SecondFinDiffPlan)
        maxerr = 5e-5

        @test all((ix, iy) for ix in 1:nx, iy in 1:ny) do (ix, iy)
            ∂x, ∂y = differentiate(vals, ix, iy, (nx, ny))
            x ,  y = (ix - 1) / (nx - 1), (iy - 1) / (ny - 1)
            tx, ty = fx(x, y), fy(x, y)
            return abs(tx - ∂x) < maxerr && abs(ty - ∂y) < maxerr
        end
    end
end

@safetestset "Simpson weights" begin
    using ..SecondFinDiffPlans: getpoints, SecondFinDiffPlan, getsimpweight

    @test [getsimpweight(Float64, ix, iy, (3, 5)) * 9 * 4 * 2 for iy in 1:5, ix in 1:3] ≈ [
        1  4 1;
        4 16 4;
        2  8 2;
        4 16 4;
        1  4 1]

    @test [getsimpweight(Float64, ix, iy, (5, 7)) for iy in 1:7, ix in 1:5] ≈ [
        1  4  2  4  1;
        4 16  8 16  4;
        2  8  4  8  2;
        4 16  8 16  4;
        2  8  4  8  2;
        4 16  8 16  4;
        1  4  2  4  1
    ] ./ (9 * 4 * 6)

    let f(x, y) = 1 + 2x + 3y + 4x^2 + 5x*y + 6y^2
        nx, ny = (3, 3)
        vals = getpoints(f, Float64, (nx, ny), SecondFinDiffPlan)
        weights = [getsimpweight(Float64, ix, iy, (nx, ny)) for iy in 1:ny, ix in 1:nx] |> vec
        @test sum(weights .* vals) ≈ 97 / 12 atol=10eps()
    end
end

@safetestset "compute!" begin
    using ..SecondFinDiffPlans: getpoints, SecondFinDiffPlan, getsimpweight
    using PoincareInvariants
    using LinearAlgebra: dot

    f(x, y) = [
         1 +  2*x +  3*y +  4*x^2 +  5*x*y +  6*y^2,
         7 +  8*x +  9*y + 10*x^2 + 11*x*y + 12*y^2,
        13 + 14*x + 15*y + 16*x^2 + 17*x*y + 18*y^2,
        19 + 20*x + 21*y + 22*x^2 + 23*x*y + 24*y^2]

    fx(x, y) = [
         2 +  8*x +  5*y,
         8 + 20*x + 11*y,
        14 + 32*x + 17*y,
        20 + 44*x + 23*y]

    fy(x, y) = [
         3 +  5*x + 12*y,
         9 + 11*x + 24*y,
        15 + 17*x + 36*y,
        21 + 23*x + 48*y]

    ω(z, ::Any, ::Any) = [0 0 z[1] z[2]; 0 0 z[3] z[4]; -z[1] -z[3] 0 0; -z[2] -z[4] 0 0]
    nx, ny = 11, 17
    integrand = getpoints(Float64, (nx, ny), SecondFinDiffPlan) do x, y
        dot(fy(x, y), ω(f(x, y), 0, nothing), fx(x, y))
    end

    ws = [getsimpweight(Float64, ix, iy, (nx, ny)) for iy in 1:ny, ix in 1:nx] |> vec
    testI = dot(ws, integrand)

    pinv = SecondPoincareInvariant{Float64, 4}(ω, (nx, ny), SecondFinDiffPlan)

    @test pinv.plan isa SecondFinDiffPlan

    points = getpoints(f, Float64, (nx, ny), SecondFinDiffPlan)
    I = compute!(pinv, points, 0, nothing)

    @test I ≈ testI atol=5e-11
end
