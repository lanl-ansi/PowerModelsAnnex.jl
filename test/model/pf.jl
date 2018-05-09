
function run_ac_pf_model(data, solver)
    model = post_ac_pf(data, Model(solver=solver))
    status = solve(model)
    return status, model
end

@testset "test ac polar pf" begin
    @testset "case $(case_file)" for (name, case_file) in case_files
        data = PMs.parse_file(case_file)
        pf_status, pf_model = run_ac_pf_model(data, ipopt_solver)
        pm_result = run_ac_pf(data, ipopt_solver)
        pm_sol = pm_result["solution"]

        base_mva = data["baseMVA"]

        for (i, bus) in data["bus"]
            if bus["bus_type"] != 4
                index = parse(Int, i)
                #println("$i, $(getvalue(pf_model[:va][index])), $(pm_sol["bus"][i]["va"])")
                @test isapprox(getvalue(pf_model[:va][index]), pm_sol["bus"][i]["va"]; atol = 1e-8)
                @test isapprox(getvalue(pf_model[:vm][index]), pm_sol["bus"][i]["vm"])
            end
        end

        for (i, gen) in data["gen"]
            if gen["gen_status"] != 0
                index = parse(Int, i)
                #println("$i, $(getvalue(pf_model[:pg][index])), $(pm_sol["gen"][i]["pg"])")
                @test isapprox(getvalue(pf_model[:pg][index]), pm_sol["gen"][i]["pg"]; atol = 1e-8)
                # multiple generators at one bus can cause this to be non-unqiue
                #@test isapprox(getvalue(pf_model[:qg][index]), pm_sol["gen"][i]["qg"])
            end
        end
    end
end

function run_soc_pf_model(data, solver)
    model = post_soc_pf(data, Model(solver=solver))
    status = solve(model)
    return status, model
end

@testset "test soc w pf" begin
    @testset "case $(case_file)" for (name, case_file) in case_files
        data = PMs.parse_file(case_file)
        pf_status, pf_model = run_soc_pf_model(data, ipopt_solver)
        pm_result = run_pf(data, SOCWRPowerModel, ipopt_solver)
        pm_sol = pm_result["solution"]

        #println(pf_status)
        #println(pm_result["status"])

        base_mva = data["baseMVA"]

        for (i, bus) in data["bus"]
            if bus["bus_type"] != 4
                index = parse(Int, i)
                #println("$i, $(getvalue(pf_model[:va][index])), $(pm_sol["bus"][i]["va"])")
                #@test isapprox(getvalue(pf_model[:va][index]), pm_sol["bus"][i]["va"]; atol = 1e-8)
                #println("$i, $(getvalue(pf_model[:w][index])), $(pm_sol["bus"][i]["vm"]^2)")
                @test isapprox(getvalue(pf_model[:w][index]), pm_sol["bus"][i]["vm"]^2; atol = 1e-3)
            end
        end

        for (i, gen) in data["gen"]
            if gen["gen_status"] != 0
                index = parse(Int, i)
                #println("$i, $(getvalue(pf_model[:pg][index])), $(pm_sol["gen"][i]["pg"])")
                @test isapprox(getvalue(pf_model[:pg][index]), pm_sol["gen"][i]["pg"]; atol = 1e-1)
                # multiple generators at one bus can cause this to be non-unqiue
                #@test isapprox(getvalue(pf_model[:qg][index]), pm_sol["gen"][i]["qg"])
            end
        end
    end
end

function run_dc_pf_model(data, solver)
    model = post_dc_pf(data, Model(solver=solver))
    status = solve(model)
    return status, model
end

@testset "test dc polar pf" begin
    @testset "case $(case_file)" for (name, case_file) in case_files
        data = PMs.parse_file(case_file)
        pf_status, pf_model = run_dc_pf_model(data, ipopt_solver)
        pm_result = run_dc_pf(data, ipopt_solver)
        pm_sol = pm_result["solution"]

        #println(pf_status)
        #println(pm_result["status"])

        # needed becouse some test networks are not DC feasible
        if pm_result["status"] == :LocalOptimal
            @test pf_status == :Optimal

            base_mva = data["baseMVA"]

            for (i, bus) in data["bus"]
                if bus["bus_type"] != 4
                    index = parse(Int, i)
                    #println("$i, $(getvalue(pf_model[:va][index])), $(pm_sol["bus"][i]["va"])")
                    @test isapprox(getvalue(pf_model[:va][index]), pm_sol["bus"][i]["va"]; atol = 1e-8)
                end
            end

            for (i, gen) in data["gen"]
                if gen["gen_status"] != 0
                    index = parse(Int, i)
                    #println("$i, $(getvalue(pf_model[:pg][index])), $(pm_sol["gen"][i]["pg"])")
                    @test isapprox(getvalue(pf_model[:pg][index]), pm_sol["gen"][i]["pg"])
                end
            end
        else
            @test pf_status == :Infeasible
            @test pm_result["status"] == :LocalInfeasible
        end
    end
end

