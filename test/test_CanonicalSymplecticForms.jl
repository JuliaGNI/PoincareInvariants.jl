@safetestset "canonical_one_form" begin
    using PoincareInvariants: canonical_one_form
    using LinearAlgebra: dot
    using Random: rand

    @test_throws ArgumentError canonical_one_form([1, 2, 3], 0.5, nothing)
    @test_throws ArgumentError canonical_one_form([0.3], 3.4, nothing)

    for mid in [3, 9, 22]
        n = mid * 2

        q = rand(mid)
        p = rand(mid)
        z = [q; p]

        @inferred canonical_one_form(z, 0.1, nothing)

        θ = canonical_one_form(z, 1.3, nothing)

        @test θ[1:mid] == p
        @test θ[mid+1:end] == zeros(mid)

        @test dot(θ, z) ≈ dot(z, θ) ≈ dot(p, q)

        v = rand(n)
        @test dot(θ, v) ≈ dot(v, θ) ≈ dot(p, v[1:mid])
    end
end

@safetestset "canonical_two_form" begin
    using PoincareInvariants: CanonicalSymplecticMatrix, canonical_two_form
    using LinearAlgebra: dot
    using Random: rand

    @test CanonicalSymplecticMatrix{Int}(2) == [0 -1; 1 0]
    @test CanonicalSymplecticMatrix{Int}(4) == [0 0 -1 0; 0 0 0 -1; 1 0 0 0; 0 1 0 0]
    @test CanonicalSymplecticMatrix{Int}(6) == [0  0  0 -1  0  0;
                                                0  0  0  0 -1  0;
                                                0  0  0  0  0 -1;
                                                1  0  0  0  0  0;
                                                0  1  0  0  0  0;
                                                0  0  1  0  0  0]

    @test_throws ArgumentError CanonicalSymplecticMatrix{Float64}(-2)
    @test_throws ArgumentError CanonicalSymplecticMatrix{Int}(0)

    function test_canonicalmatrix(mat, n, ::Type{T}) where T
        mid = n ÷ 2

        for i in 1:n, j in 1:n
            if i == j - mid
                mat[i, j] === T(-1) || return false
            elseif i - mid == j
                mat[i, j] === T(1) || return false
            else
                mat[i, j] === T(0) || return false
            end
        end

        return true
    end

    @testset "CanonicalSymplecticMatrix{$T}($n)" for T in [Int, Float64], n in [10, 64, 1002]
        @test CanonicalSymplecticMatrix(n) === CanonicalSymplecticMatrix{Int}(n)

        C = CanonicalSymplecticMatrix{T}(n)
        @test size(C) == (n, n)
        @test eltype(C) == T

        mid = n ÷ 2

        @test test_canonicalmatrix(C, n, T)

        @test_throws BoundsError C[0, 4]
        @test_throws BoundsError C[8, -2]
        @test_throws BoundsError C[2    , n + 1]
        @test_throws BoundsError C[n + 1,     3]

        @test C * collect(1:n) == [(-mid-1:-1:-n)..., 1:mid...]

        @test dot(ones(n), C, ones(n)) == 0

        v = rand(n); w = rand(n)
        vCw = dot(v, C, w)

        @test vCw ≈ dot(v[mid+1:end], w[1:mid]) - dot(v[1:mid], w[mid+1:end])
        @test vCw ≈ dot(v, C * w)

        @test_throws ArgumentError CanonicalSymplecticMatrix{T}(n + 1)

        @inferred CanonicalSymplecticMatrix{Float64} canonical_two_form(v, 0.1, nothing)

        @test canonical_two_form(v, 0.1, nothing) == CanonicalSymplecticMatrix(n)
    end
end
