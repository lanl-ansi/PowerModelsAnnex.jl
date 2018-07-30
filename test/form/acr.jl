
@testset "test acr nl" begin
    @testset "3-bus case" begin
        pm = build_generic_model(case_files["case3"], NLACRPowerModel, PowerModels.post_opf)

        result = solve_generic_model(pm, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 5907; atol = 1e0)
    end
    @testset "5-bus pjm case" begin
        pm = build_generic_model(case_files["case5"], NLACRPowerModel, PowerModels.post_opf)

        result = solve_generic_model(pm, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 18269; atol = 1e0) # Increasing tolerance
        # here as the result seems to depend on the machine
    end
    @testset "30-bus ieee case" begin
        pm = build_generic_model(case_files["case30"], NLACRPowerModel, PowerModels.post_opf)

        result = solve_generic_model(pm, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 204.9; atol = 1e0)
    end
end
