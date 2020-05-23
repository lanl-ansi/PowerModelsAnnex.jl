
@testset "test ac api" begin
    @testset "3-bus case" begin
        result = run_opf_api(case_files["case3"], ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 1.3371; atol = 1e-3)
        @test isapprox(result["solution"]["load"]["1"]["pd"], 1.471; atol = 1e-2)
        @test isapprox(result["solution"]["load"]["1"]["qd"], 0.4; atol = 1e-2)
    end
    # started failing 05/22/2020 when ipopt moved to jll artifacts
    # @testset "5-bus pjm case" begin
    #     result = run_opf_api(case_files["case5"], ACPPowerModel, ipopt_solver)

    #     @test result["termination_status"] == LOCALLY_SOLVED
    #     @test isapprox(result["objective"], 2.6872; atol = 1e-3)
    #     @test isapprox(result["solution"]["load"]["3"]["pd"], 10.754; atol = 1e-2)
    #     @test isapprox(result["solution"]["load"]["3"]["qd"], 1.3147; atol = 1e-2)
    # end
    @testset "14-bus ieee case" begin
        result = run_opf_api(case_files["case14"], ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 3.98769; atol = 1e-3)
        @test isapprox(result["solution"]["load"]["1"]["pd"], 0.8653; atol = 1e-2)
        @test isapprox(result["solution"]["load"]["1"]["qd"], 0.127; atol = 1e-2)
    end
    @testset "30-bus ieee case" begin
        result = run_opf_api(case_files["case30"], ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 1.6628; atol = 1e-3)
        @test isapprox(result["solution"]["load"]["1"]["pd"], 0.361; atol = 1e-2)
        @test isapprox(result["solution"]["load"]["1"]["qd"], 0.127; atol = 1e-2)
    end
end
