using MultiObjectiveProblems
using Test

@testset "MultiObjectiveProblems.jl" begin
    T1 = MHT1();
    @test num_objectives(T1) == num_vars(T1) == 2
    ps = get_pareto_set( T1 );
    ps_p1 = get_points( ps, 10 );
    ps_p2 = get_points( ps, 10; method = :random)
    @test length(ps_p1) == length(ps_p2)
    @test all( [sum(p) ≈ 10 for p ∈ ps_p1] )
    @test all( [sum(p) ≈ 10 for p ∈ ps_p2] )
    @test get_starting_point(T1,1) == [-34.24; 47.06]
end