@safetestset "CanonicalSymplecticTwoForm" begin
    using PoincareInvariants: CanonicalSymplecticTwoForm
    using LinearAlgebra: dot
    using Random: rand

    @test CanonicalSymplecticTwoForm{Int}(2) == [0 -1; 1 0]
    @test CanonicalSymplecticTwoForm{Int}(4) == [0 0 -1 0; 0 0 0 -1; 1 0 0 0; 0 1 0 0]
    @test CanonicalSymplecticTwoForm{Int}(6) == [0  0  0 -1  0  0;
                                                0  0  0  0 -1  0;
                                                0  0  0  0  0 -1;
                                                1  0  0  0  0  0;
                                                0  1  0  0  0  0;
                                                0  0  1  0  0  0]

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

    @testset "CanonicalSymplecticTwoForm{$T}($n)" for T in [Int, Float64], n in [10, 64, 1002]
        @test CanonicalSymplecticTwoForm(n) === CanonicalSymplecticTwoForm{Int}(n)

        C = CanonicalSymplecticTwoForm{T}(n)
        @test size(C) == (n, n)
        @test eltype(C) == T

        mid = n ÷ 2

        @test test_canonicalmatrix(C, n, T)

        @test C * collect(1:n) == [(-mid-1:-1:-n)..., 1:mid...]

        @test dot(ones(n), C, ones(n)) == 0

        v = rand(n); w = rand(n)
        vCw = dot(v, C, w)

        @test vCw ≈ dot(v[mid+1:end], w[1:mid]) - dot(v[1:mid], w[mid+1:end])
        @test vCw ≈ dot(v, C * w)

        @test_throws ArgumentError CanonicalSymplecticTwoForm{T}(n + 1)
    end
end
