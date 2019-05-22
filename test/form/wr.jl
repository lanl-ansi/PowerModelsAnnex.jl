
@testset "test wr oa" begin
    @testset "3-bus case" begin
        pm = PMs.build_generic_model(case_files["case3"], SOCWROAPowerModel, PowerModels.post_opf)

        p = PMs.var(pm, :p)
        for v in values(p)
            JuMP.set_start_value(v, 1.0)
        end
        q = PMs.var(pm, :q)
        for v in values(q)
            JuMP.set_start_value(v, 1.0)
        end

        result = PMs.solve_generic_model(pm, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 5746.7; atol = 1e0)
    end
    @testset "5-bus pjm case" begin
        pm = PMs.build_generic_model(case_files["case5"], SOCWROAPowerModel, PowerModels.post_opf)

        p = PMs.var(pm, :p)
        for v in values(p)
            JuMP.set_start_value(v, 1.0)
        end
        q = PMs.var(pm, :q)
        for v in values(q)
            JuMP.set_start_value(v, 1.0)
        end

        result = PMs.solve_generic_model(pm, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 14999.71; atol = 1e2) # Increasing tolerance
        # here as the result seems to depend on the machine
    end
    @testset "30-bus ieee case" begin
        pm = PMs.build_generic_model(case_files["case30"], SOCWROAPowerModel, PowerModels.post_opf)

        p = PMs.var(pm, :p)
        for v in values(p)
            JuMP.set_start_value(v, 1.0)
        end
        q = PMs.var(pm, :q)
        for v in values(q)
            JuMP.set_start_value(v, 1.0)
        end

        result = PMs.solve_generic_model(pm, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 172.41; atol = 1e0)
    end
end



@testset "test qc tri without linking" begin
    @testset "3-bus case" begin
        pm = PMs.build_generic_model(case_files["case3"], QCWRTriNoLinkPowerModel, PowerModels.post_opf)

        result = PMs.solve_generic_model(pm, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 5817.58; atol = 1e0)
    end
    @testset "5-bus pjm case" begin
        pm = PMs.build_generic_model(case_files["case5"], QCWRTriNoLinkPowerModel, PowerModels.post_opf)

        result = PMs.solve_generic_model(pm, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 15051.6; atol = 1e2)
    end
    @testset "30-bus ieee case" begin
        pm = PMs.build_generic_model(case_files["case30"], QCWRTriNoLinkPowerModel, PowerModels.post_opf)

        result = PMs.solve_generic_model(pm, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 173.806; atol = 1e0)
    end
end
