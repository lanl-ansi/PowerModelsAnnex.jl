
@testset "test wr oa" begin
    @testset "3-bus case" begin
        pm = build_generic_model(case_files["case3"], SOCWROAPowerModel, PowerModels.post_opf)

        p = var(pm, :p)
        for k in keys(p)
            setvalue(p[k[1]], 1.0)
        end
        q = var(pm, :q)
        for k in keys(q)
            setvalue(q[k[1]], 1.0)
        end

        result = solve_generic_model(pm, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 5746.7; atol = 1e0)
    end
    @testset "5-bus pjm case" begin
        pm = build_generic_model(case_files["case5"], SOCWROAPowerModel, PowerModels.post_opf)

        p = var(pm, :p)
        for k in keys(p)
            setvalue(p[k[1]], 1.0)
        end
        q = var(pm, :q)
        for k in keys(q)
            setvalue(q[k[1]], 1.0)
        end

        result = solve_generic_model(pm, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 14999.71; atol = 1e2) # Increasing tolerance
        # here as the result seems to depend on the machine
    end
    @testset "30-bus ieee case" begin
        pm = build_generic_model(case_files["case30"], SOCWROAPowerModel, PowerModels.post_opf)

        p = var(pm, :p)
        for k in keys(p)
            setvalue(p[k[1]], 1.0)
        end
        q = var(pm, :q)
        for k in keys(q)
            setvalue(q[k[1]], 1.0)
        end

        result = solve_generic_model(pm, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 172.41; atol = 1e0)
    end
end



@testset "test qc tri without linking" begin
    @testset "3-bus case" begin
        pm = build_generic_model(case_files["case3"], QCWRTriNoLinkPowerModel, PowerModels.post_opf)

        result = solve_generic_model(pm, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 5817.58; atol = 1e0)
    end
    @testset "5-bus pjm case" begin
        pm = build_generic_model(case_files["case5"], QCWRTriNoLinkPowerModel, PowerModels.post_opf)

        result = solve_generic_model(pm, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 15051.6; atol = 1e2) 
    end
    @testset "30-bus ieee case" begin
        pm = build_generic_model(case_files["case30"], QCWRTriNoLinkPowerModel, PowerModels.post_opf)

        result = solve_generic_model(pm, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 173.806; atol = 1e0)
    end
end
