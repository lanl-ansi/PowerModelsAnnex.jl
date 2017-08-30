
@testset "test ac sad" begin
    @testset "3-bus case" begin
        result = PowerModelsAnnex.run_sad_opf("../../PowerModels/test/data/case3.m", ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.3114; atol = 1e-2)
    end
    @testset "5-bus pjm case" begin
        result = PowerModelsAnnex.run_sad_opf("../../PowerModels/test/data/case5.m", ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.02211; atol = 1e-2)
    end
    @testset "30-bus ieee case" begin
        result = PowerModelsAnnex.run_sad_opf("../../PowerModels/test/data/case30.m", ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.1522; atol = 1e-2)
    end
end
