@safetestset "CanonicalSymplecticMatrix" begin
    using PoincareInvariants: CanonicalSymplecticMatrix
    using LinearAlgebra: dot

    @test CanonicalSymplecticMatrix{Int}(2) == [0 -1; 1 0]
    @test CanonicalSymplecticMatrix{Int}(4) == [0 0 -1 0; 0 0 0 -1; 1 0 0 0; 0 1 0 0]
    @test CanonicalSymplecticMatrix{Int}(6) == [0  0  0 -1  0  0;
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

    @testset "CanonicalSymplecticMatrix{$T}($n)" for T in [Int, Float64], n in [10, 100, 1000]
        @test CanonicalSymplecticMatrix(n) === CanonicalSymplecticMatrix{Int}(n)

        C = CanonicalSymplecticMatrix{T}(n)
        @test size(C) == (n, n)
        @test eltype(C) == T

        mid = n ÷ 2

        @test test_canonicalmatrix(C, n, T)

        @test C * collect(1:n) == [(-mid-1:-1:-n)..., 1:mid...]

        @test dot(ones(n), C, ones(n)) == 0

        @test_throws ArgumentError CanonicalSymplecticMatrix{T}(n + 1)
    end
end
