@safetestset "SecondPoincareInvariant Constructors" begin
    using PoincareInvariants

    for D in [2, 10], T in [Float32, Float64], N in [1, 2, 5, 10, 100]
        # SecondPoincareInvariant{T}(Ω::AbstractMatrix, D::Integer, N::Integer)
        Ω = rand(Int, D, D)

        let pinv = SecondPoincareInvariant{D, T}(Ω, N)
            @test pinv isa SecondPoincareInvariant{D, T, Matrix{Int}, <:Any}
            @test pinv.Ω === Ω
        end

        @test_throws ArgumentError SecondPoincareInvariant{D + 1, T}(rand(Int, D - 1, D), N)

        @test_throws ArgumentError SecondPoincareInvariant{D, T}(rand(Int, D, D), 0)
        @test_throws ArgumentError SecondPoincareInvariant{D, T}(rand(Int, D, D), -N)

        # SecondPoincareInvariant{D, T}(Ω::Function, N::Integer, Val(true))
        ω!(out, z, p, t) = copyto!(out, Ω)
        
        let pinv = SecondPoincareInvariant{D, T}(ω!, N)
            @test pinv isa SecondPoincareInvariant{D, T, typeof(ω!), <:Any}
            @test pinv.Ω === ω!
        end

        let pinv = SecondPoincareInvariant{D, T}(ω!, N, Val(true))
            @test pinv isa SecondPoincareInvariant{D, T, typeof(ω!), <:Any}
            @test pinv.Ω === ω!
        end

        # SecondPoincareInvariant{D, T}(Ω::Function, N::Integer, Val(false))
        ω(z, p, t) = copy(Ω)
        
        let pinv = SecondPoincareInvariant{D, T}(ω, N)
            @test pinv isa SecondPoincareInvariant{D, T, typeof(ω), <:Any}
            @test pinv.Ω === ω
        end

        let pinv = SecondPoincareInvariant{D, T}(ω, N, Val(false))
            @test pinv isa SecondPoincareInvariant{D, T, typeof(ω), <:Any}
            @test pinv.Ω === ω
        end
    end
end

@safetestset "getpoints and getpointnum" begin
    using PoincareInvariants

    for N in [1, 2, 5, 10, 100], T in [Float32, Float64]
        pinv = SecondPoincareInvariant{6, T}(rand(Int, 6, 6), N)
        points = getpoints(pinv)

        @test points isa Matrix{T}
        @test N ≤ size(points, 1) ≤ 2N
        @test getpointnum(pinv) == size(points, 1)
        @test size(points, 2) == 2

        let ω(z, p, t) = rand(Int, 2, 2)
            pinv = SecondPoincareInvariant{2, T}(ω, N, Val(false))
            
            @test getpoints(pinv) == points
            @test getpointnum(pinv) == size(points, 1)
        end

        let ω!(out, z, p, t) = copyto!(out, rand(Int, 4, 4))
            pinv = SecondPoincareInvariant{4, T}(ω!, N, Val(true))
            
            @test getpoints(pinv) == points
            @test getpointnum(pinv) == size(points, 1)
        end
    end
end

@safetestset "compute" begin
    using PoincareInvariants

    @testset "Constant Ω" for D in [2, 10], T in [Float32, Float64], N in [10, 100]
        Ω = CanonicalSymplecticMatrix(D)
        pinv = SecondPoincareInvariant{D, T}(Ω, N)

        parampoints = getpoints(pinv)

        phasepoints = zeros(T, getpointnum(pinv), D)
        phasepoints[:, 1] = parampoints[:, 1]
        phasepoints[:, N ÷ 2 + 1] = parampoints[:, 2]

        out = compute(pinv, phasepoints)
        @test out isa T
        @test out ≈ 1 rtol=√eps(T)
    end

    @testset "ω(z, (a = 5, b = 3.0), 9.5)" for D in [2, 10], T in [Float32, Float64], N in [10, 100]
        function ω(z, p, t)
            p === (a = 5, b = 3.0) || error()
            t === 9.5 || error()
            CanonicalSymplecticMatrix(D)
        end

        pinv = SecondPoincareInvariant{D, T}(ω, N)

        parampoints = getpoints(pinv)

        phasepoints = zeros(T, getpointnum(pinv), D)
        phasepoints[:, 1] = parampoints[:, 1]
        phasepoints[:, N ÷ 2 + 1] = parampoints[:, 2]

        out = compute(pinv, phasepoints, (a = 5, b = 3.0), 9.5)
        @test out isa T
        @test out ≈ 1 rtol=√eps(T)
    end

    @testset "ω!(out, z, 39, 1f0)" for D in [2, 10], T in [Float32, Float64], N in [10, 100]
        function ω!(out, z, p, t)
            p === 39 || error()
            t === 1f0 || error()
            cpoyto!(out, CanonicalSymplecticMatrix(D))
        end

        pinv = SecondPoincareInvariant{D, T}(ω!, N)

        parampoints = getpoints(pinv)

        phasepoints = zeros(T, getpointnum(pinv), D)
        phasepoints[:, 1] = parampoints[:, 1]
        phasepoints[:, N ÷ 2 + 1] = parampoints[:, 2]

        out = compute(pinv, phasepoints)
        @test out isa T
        @test out ≈ 1 rtol=√eps(T)
    end
end
