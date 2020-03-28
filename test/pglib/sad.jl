
@testset "test ac sad" begin
    @testset "3-bus case" begin
        result = PowerModelsAnnex.run_opf_sad(case_files["case3"], PMs.ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.3114; atol = 1e-2)
    end
    @testset "5-bus pjm case" begin
        result = PowerModelsAnnex.run_opf_sad(case_files["case5"], PMs.ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.02211; atol = 1e-2)
    end
    @testset "14-bus ieee case" begin
        result = PowerModelsAnnex.run_opf_sad(case_files["case14"], PMs.ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.033097; atol = 1e-2)
    end
    @testset "30-bus ieee case" begin
        result = PowerModelsAnnex.run_opf_sad(case_files["case30"], PMs.ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.1522; atol = 1e-2)
    end
end
