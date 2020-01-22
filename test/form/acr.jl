
@testset "test acr nl" begin
    @testset "3-bus case" begin
        pm = PMs.instantiate_model(case_files["case3"], NLACRPowerModel, PMs.build_opf)

        result = PMs.optimize_model!(pm, optimizer=ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 5907; atol = 1e0)
    end
    @testset "5-bus pjm case" begin
        pm = PMs.instantiate_model(case_files["case5"], NLACRPowerModel, PMs.build_opf)

        result = PMs.optimize_model!(pm, optimizer=ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 18269; atol = 1e0) # Increasing tolerance
        # here as the result seems to depend on the machine
    end
    @testset "30-bus ieee case" begin
        pm = PMs.instantiate_model(case_files["case30"], NLACRPowerModel, PMs.build_opf)

        result = PMs.optimize_model!(pm, optimizer=ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 204.9; atol = 1e0)
    end
end
