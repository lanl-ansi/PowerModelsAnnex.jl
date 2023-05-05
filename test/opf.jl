@testset "test ac cop" begin
    @testset "3-bus case" begin
        result = solve_opf_cop(case_files["case3"], ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol = 1e-3)
        @test isapprox(result["solution"]["gen"]["1"]["pg"], 1.581; atol = 1e-2)
        @test isapprox(result["solution"]["gen"]["2"]["pg"], 1.599; atol = 1e-2)
        @test isapprox(result["solution"]["bus"]["1"]["vm"], 1.099; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["vm"], 0.926; atol = 1e-3)
    end
    @testset "5-bus pjm case" begin
        result = solve_opf_cop(case_files["case5"], ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.136996; atol = 1e-3)
        @test isapprox(result["solution"]["gen"]["1"]["pg"], 0.316; atol = 1e-2)
        @test isapprox(result["solution"]["gen"]["2"]["pg"], 1.616; atol = 1e-2)
        @test isapprox(result["solution"]["bus"]["1"]["vm"], 0.977; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["vm"], 0.975; atol = 1e-3)
    end
    @testset "14-bus ieee case" begin
        result = solve_opf_cop(case_files["case14"], ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0024798; atol = 1e-3)
        @test isapprox(result["solution"]["gen"]["1"]["pg"], 2.323; atol = 1e-2)
        @test isapprox(result["solution"]["gen"]["2"]["pg"], 0.399; atol = 1e-2)
        @test isapprox(result["solution"]["bus"]["1"]["vm"], 1.059; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["vm"], 1.037; atol = 1e-3)
    end
    @testset "30-bus ieee case" begin
        result = solve_opf_cop(case_files["case30"], ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol = 1e-3)
        @test isapprox(result["solution"]["gen"]["1"]["pg"], 2.188; atol = 1e-2)
        @test isapprox(result["solution"]["gen"]["2"]["pg"], 0.800; atol = 1e-2)
        @test isapprox(result["solution"]["bus"]["1"]["vm"], 1.059; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["vm"], 1.036; atol = 1e-3)
    end
end

@testset "test dc cop" begin
    @testset "5-bus pjm case" begin
        result = solve_opf_cop(case_files["case5"], DCPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.00717; atol = 1e-3)
        @test isapprox(result["solution"]["gen"]["1"]["pg"], 0.369; atol = 1e-2)
        @test isapprox(result["solution"]["gen"]["2"]["pg"], 1.669; atol = 1e-2)
    end
end


