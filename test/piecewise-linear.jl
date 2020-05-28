@testset "test ac polar pwl opf" begin
    @testset "5-bus with pwl delta model" begin
        result = run_opf_pwl_delta(case_file_pwl, ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 42895; atol = 1e0)
    end
    @testset "5-bus with pwl lambda model" begin
        result = run_opf_pwl_lambda(case_file_pwl, ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 42895; atol = 1e0)
    end
    @testset "5-bus with pwl phi model" begin
        result = run_opf_pwl_phi(case_file_pwl, ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 42895; atol = 1e0)
    end
    @testset "5-bus with pwl psi model" begin
        result = run_opf_pwl_psi(case_file_pwl, ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 42895; atol = 1e0)
    end
end

@testset "test soc pwl opf" begin
    @testset "5-bus with pwl delta model" begin
        result = run_opf_pwl_delta(case_file_pwl, SOCWRPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 42895; atol = 1e0)
    end
    @testset "5-bus with pwl lambda model" begin
        result = run_opf_pwl_lambda(case_file_pwl, SOCWRPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 42895; atol = 1e0)
    end
    @testset "5-bus with pwl phi model" begin
        result = run_opf_pwl_phi(case_file_pwl, SOCWRPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 42895; atol = 1e0)
    end
    @testset "5-bus with pwl psi model" begin
        result = run_opf_pwl_psi(case_file_pwl, SOCWRPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 42895; atol = 1e0)
    end
end

@testset "test dc pwl opf" begin
    @testset "5-bus with pwl delta model" begin
        result = run_opf_pwl_delta(case_file_pwl, DCPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 42565; atol = 1e0)
    end
    @testset "5-bus with pwl lambda model" begin
        result = run_opf_pwl_lambda(case_file_pwl, DCPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 42565; atol = 1e0)
    end
    @testset "5-bus with pwl phi model" begin
        result = run_opf_pwl_phi(case_file_pwl, DCPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 42565; atol = 1e0)
    end
    @testset "5-bus with pwl psi model" begin
        result = run_opf_pwl_psi(case_file_pwl, DCPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 42565; atol = 1e0)
    end
end

