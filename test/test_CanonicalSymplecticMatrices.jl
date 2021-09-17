using PoincareInvariants.CanonicalSymplecticMatrices
using LinearAlgebra: dot

@test CanonicalSymplecticMatrix{Int}(2) == [0 -1; 1 0]
@test CanonicalSymplecticMatrix{Int}(4) == [0 0 -1 0; 0 0 0 -1; 1 0 0 0; 0 1 0 0]
@test CanonicalSymplecticMatrix{Int}(6) == [0 0 0 -1 0 0;
                                            0 0 0 0 -1 0;
                                            0 0 0 0 0 -1;
                                            1 0 0 0 0 0;
                                            0 1 0 0 0 0;
                                            0 0 1 0 0 0]

@testset "CanonicalSymplecticMatrix{$T}($n)" for T in [Int, Float64], n in [2:2:10..., 20:10:100...]
    @test CanonicalSymplecticMatrix(n) === CanonicalSymplecticMatrix{Int}(n)

    C = CanonicalSymplecticMatrix{T}(n)
    @test size(C) == (n, n)
    @test eltype(C) == T

    mid = n ÷ 2

    for i in 1:n, j in 1:n
        if i == j - mid
            @test C[i, j] === T(-1)
        elseif i - mid == j
            @test C[i, j] === T(1)
        else
            @test C[i, j] === T(0)
        end
    end

    @test C * collect(1:n) == [(-mid-1:-1:-n)..., 1:mid...]

    @test dot(ones(n), C, ones(n)) == 0

    @test_throws ArgumentError CanonicalSymplecticMatrix{T}(n + 1)
end
