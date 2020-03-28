
@testset "test ac sad" begin
    @testset "3-bus case" begin
        result = run_opf_sad(case_files["case3"], ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.3114; atol = 1e-2)
    end
    @testset "5-bus pjm case" begin
        result = run_opf_sad(case_files["case5"], ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.02211; atol = 1e-2)
    end
    @testset "14-bus ieee case" begin
        result = run_opf_sad(case_files["case14"], ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.033097; atol = 1e-2)
    end
    @testset "30-bus ieee case" begin
        result = run_opf_sad(case_files["case30"], ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.1522; atol = 1e-2)
    end
end
