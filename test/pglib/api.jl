
@testset "test ac api" begin
    @testset "3-bus case" begin
        result = PowerModelsAnnex.run_api_opf("../../PowerModels/test/data/matpower/case3.m", ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 1.3371; atol = 1e-3)
        @test isapprox(result["solution"]["load"]["1"]["pd"], 1.471; atol = 1e-2)
    end
    @testset "5-bus pjm case" begin
        result = PowerModelsAnnex.run_api_opf("../../PowerModels/test/data/matpower/case5.m", ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 2.6872; atol = 1e-3)
        @test isapprox(result["solution"]["load"]["3"]["pd"], 10.754; atol = 1e-2)
    end
    @testset "30-bus ieee case" begin
        result = PowerModelsAnnex.run_api_opf("../../PowerModels/test/data/matpower/case30.m", ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 1.6628; atol = 1e-3)
        @test isapprox(result["solution"]["load"]["1"]["pd"], 0.361; atol = 1e-2)
    end
end
