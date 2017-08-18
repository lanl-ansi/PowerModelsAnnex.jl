
function run_ac_opf_model(data, solver)
    model = post_ac_opf(data, Model(solver=solver))
    solve(model)
    return model
end

@testset "test ac polar opf" begin
    @testset "case $(case_file)" for case_file in case_files
        data = PMs.parse_file(case_file)
        opf_model = run_ac_opf_model(data, ipopt_solver)
        pm_result = run_ac_opf(data, ipopt_solver)
        pm_sol = pm_result["solution"]

        base_mva = data["baseMVA"]

        for (i, bus) in data["bus"]
            if bus["bus_type"] != 4
                index = parse(Int, i)
                #println("$i, $(getvalue(opf_model[:t][index])), $(pm_sol["bus"][i]["va"]*pi/180)")
                @test isapprox(getvalue(opf_model[:t][index]), pm_sol["bus"][i]["va"]*pi/180; atol = 1e-8)
                @test isapprox(getvalue(opf_model[:v][index]), pm_sol["bus"][i]["vm"])
            end
        end

        for (i, gen) in data["gen"]
            if gen["gen_status"] != 0
                index = parse(Int, i)
                @test isapprox(getvalue(opf_model[:pg][index]), pm_sol["gen"][i]["pg"]/base_mva)
                # multiple generators at one bus can cause this to be non-unqiue
                #@test isapprox(getvalue(opf_model[:qg][index]), pm_sol["gen"][i]["qg"]/base_mva)
            end
        end
    end
end