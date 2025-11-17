using CopositivityDiscriminants
using Test
using HomotopyContinuation
using Oscar

@testset "CopositivityDiscriminants.jl" begin
    
@testset "check_copositivity — nonseparable" begin
    HomotopyContinuation.@var x[1:4]
    f = x[1]^40 + x[2]^40 + x[3]^40 + x[4]^40 -
        float((10/9)^(9/10) * 40^(1/10)) * x[1]*x[2]*x[3]*x[4] + 1

    r = check_copositivity(f, true)

    @test r.method == :nonseparable
    @test r.copositive === missing
    @test r.t_min ≈ 1.0 atol=1e-12
end
@testset "check_copositivity" begin
    HomotopyContinuation.@var x[1:4]
    f = x[1]^40 + x[2]^40 + x[3]^40 + x[4]^40 -
        float((10/9)^(9/10) * 40^(1/10)) * x[1]*x[2]*x[3]*x[4] + 1

    r = check_copositivity(f, false)

    @test r.copositive === missing
    @test r.t_min ≈ 1.0 atol=1e-12
end


end
