@safetestset "Padua" begin include("test_Padua.jl") end

# @safetestset "Coefficient and Point Counts" begin
#     using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation:
#         getdegree, getcoeffnum, getfullpointnum, getpaduanum, checkpaduanum

#     for n in [1:100..., 135, 752, 1000, 5531]
#         paduanum = getpaduanum(n)
#         @test paduanum == (n + 1) * (n + 2) ÷ 2
#         @test getdegree(paduanum) == n
#         @test (checkpaduanum(paduanum); true)
#         @test_throws ArgumentError checkpaduanum(paduanum + 1)
#         @test_throws ArgumentError checkpaduanum(paduanum - 1)
        
#         @test getcoeffnum(n) == paduanum
#         @test getfullpointnum(n - 1) == n^2
#         @test getcoeffnum(n - 1) == (getfullpointnum(n - 1) - n) ÷ 2 + n
#     end
# end

# @safetestset "getpaduapoints" begin
#     using PoincareInvariants.SecondPoincareInvariants.ChebyshevImplementation:
#         getpaduapoints, getpaduanum

#     x(m) = cospi(m / 2)
#     y(m, l, n) = iseven(m) ? cospi((2l - 1) / (n + 1)) : cospi((2l - 2) / (n + 1))

#     @testset "Points for Degree $n" for n in 2:2
#         points = getpaduapoints(n)

#         @test points isa Vector{SVector{Float64}}
#         @test length(points) == getpaduanum(n)

#         c = 1
#         for m in n:-1:0
#             lmax = fld(n, 2) + mod(m, 2) + 1
#             for l in lmax:-1:1
#                 @test points[c] ≈ [x(m), y(m, l, n)]
#                 c += 1
#             end
#         end
#     end
# end

# @safetestset "paduatransform!" begin
#     using PoincareInvariants.SecondPoincareInvariants: getpaduapoints, paduatransform!
#     using FastTransforms: plan_paduatransform!

#     f(v) = v[1] * sin(v[2])

#     paduapoints = getpaduapoints(20)
#     v = f.(eachrow(paduapoints))

#     plan = plan_paduatransform!(v, Val{false})

#     coeffs = plan * copy(v)

#     out = similar(v)
#     paduatransform!(out, v, plan)

#     @test out == coeffs
#     @test v == f.(eachrow(paduapoints))
# end
