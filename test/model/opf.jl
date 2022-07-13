
function solve_ac_opf_model(data, optimizer)
    model = build_ac_opf(data, JuMP.Model(optimizer))
    JuMP.optimize!(model)
    return JuMP.termination_status(model), model
end

@testset "test ac polar opf" begin
    @testset "case $(name)" for (name, case_file) in case_files
        data = parse_file(case_file)
        opf_status, opf_model = solve_ac_opf_model(data, ipopt_solver)
        pm_result = solve_ac_opf(data, ipopt_solver)
        pm_sol = pm_result["solution"]

        @test isapprox(JuMP.objective_value(opf_model), pm_result["objective"]; atol = 1e-5)

        base_mva = data["baseMVA"]

        for (i, bus) in data["bus"]
            if bus["bus_type"] != 4
                index = parse(Int, i)
                #println("$i, $(JuMP.value(opf_model[:va][index])), $(pm_sol["bus"][i]["va"])")
                #println("$i, $(JuMP.value(opf_model[:vm][index])), $(pm_sol["bus"][i]["vm"])")
                @test isapprox(JuMP.value(opf_model[:va][index]), pm_sol["bus"][i]["va"]; atol = 1e-8)
                @test isapprox(JuMP.value(opf_model[:vm][index]), pm_sol["bus"][i]["vm"])
            end
        end

        for (i, gen) in data["gen"]
            if gen["gen_status"] != 0
                index = parse(Int, i)
                #println("$i, $(JuMP.value(opf_model[:pg][index])), $(pm_sol["gen"][i]["pg"])")
                @test isapprox(JuMP.value(opf_model[:pg][index]), pm_sol["gen"][i]["pg"])
                # multiple generators at one bus can cause this to be non-unqiue
                #@test isapprox(JuMP.value(opf_model[:qg][index]), pm_sol["gen"][i]["qg"])
            end
        end
    end
end

function solve_soc_opf_model(data, optimizer)
    model =  build_soc_opf(data, JuMP.Model(optimizer))
    JuMP.optimize!(model)
    return JuMP.termination_status(model), model
end

@testset "test soc w opf" begin
    @testset "case $(name)" for (name, case_file) in case_files
        data = parse_file(case_file)
        opf_status, opf_model = solve_soc_opf_model(data, ipopt_solver)
        pm_result = solve_opf(data, SOCWRPowerModel, ipopt_solver)
        pm_sol = pm_result["solution"]

        @test isapprox(JuMP.objective_value(opf_model), pm_result["objective"]; atol = 1e-5)

        base_mva = data["baseMVA"]

        for (i, bus) in data["bus"]
            if bus["bus_type"] != 4
                index = parse(Int, i)
                #println("$i, $(JuMP.value(opf_model[:va][index])), $(pm_sol["bus"][i]["va"])")
                #@test isapprox(JuMP.value(opf_model[:va][index]), pm_sol["bus"][i]["va"]; atol = 1e-8)

                #println("$i, $(JuMP.value(opf_model[:w][index])), $(pm_sol["bus"][i]["vm"]^2)")
                @test isapprox(JuMP.value(opf_model[:w][index]), pm_sol["bus"][i]["w"]; atol = 1e-6)
            end
        end

        for (i, gen) in data["gen"]
            if gen["gen_status"] != 0
                index = parse(Int, i)
                #println("$i, $(JuMP.value(opf_model[:pg][index])), $(pm_sol["gen"][i]["pg"])")
                @test isapprox(JuMP.value(opf_model[:pg][index]), pm_sol["gen"][i]["pg"]; atol = 1e-6)
                # multiple generators at one bus can cause this to be non-unqiue
                #@test isapprox(JuMP.value(opf_model[:qg][index]), pm_sol["gen"][i]["qg"])
            end
        end
    end
end


function solve_qc_opf_model(data, optimizer)
    model =  build_qc_opf(data, JuMP.Model(optimizer))
    JuMP.optimize!(model)
    return JuMP.termination_status(model), model
end

@testset "test qc w+l opf" begin
    @testset "case $(name)" for (name, case_file) in case_files
        data = parse_file(case_file)
        opf_status, opf_model = solve_qc_opf_model(data, ipopt_solver)
        pm_result = solve_opf(data, QCLSPowerModel, ipopt_solver)
        pm_sol = pm_result["solution"]

        @test isapprox(JuMP.objective_value(opf_model), pm_result["objective"]; atol = 1e-5)

        base_mva = data["baseMVA"]

        for (i, bus) in data["bus"]
            if bus["bus_type"] != 4
                index = parse(Int, i)
                #println("$i, $(JuMP.value(opf_model[:va][index])), $(pm_sol["bus"][i]["va"])")
                #@test isapprox(JuMP.value(opf_model[:va][index]), pm_sol["bus"][i]["va"]; atol = 1e-8)

                #println("$i, $(JuMP.value(opf_model[:vm][index])), $(pm_sol["bus"][i]["vm"])")
                @test isapprox(JuMP.value(opf_model[:vm][index]), pm_sol["bus"][i]["vm"]; atol = 1e-5)
            end
        end

        for (i, gen) in data["gen"]
            if gen["gen_status"] != 0
                index = parse(Int, i)
                #println("$i, $(JuMP.value(opf_model[:pg][index])), $(pm_sol["gen"][i]["pg"])")
                @test isapprox(JuMP.value(opf_model[:pg][index]), pm_sol["gen"][i]["pg"]; atol = 1e-5)
                # multiple generators at one bus can cause this to be non-unqiue
                #@test isapprox(JuMP.value(opf_model[:qg][index]), pm_sol["gen"][i]["qg"])
            end
        end
    end
end


function solve_dc_opf_model(data, optimizer)
    model = build_dc_opf(data, JuMP.Model(optimizer))
    JuMP.optimize!(model)
    return JuMP.termination_status(model), model
end

@testset "test dc polar opf" begin
    @testset "case $(name)" for (name, case_file) in case_files
        data = parse_file(case_file)
        opf_status, opf_model = solve_dc_opf_model(data, ipopt_solver)
        pm_result = solve_dc_opf(data, ipopt_solver)
        pm_sol = pm_result["solution"]

        #println(opf_status)
        #println(pm_result["status"])

        @test isapprox(JuMP.objective_value(opf_model), pm_result["objective"]; atol = 1e-5)

        # needed becouse some test networks are not DC feasible
        if pm_result["termination_status"] == LOCALLY_SOLVED
            @test opf_status == LOCALLY_SOLVED

            base_mva = data["baseMVA"]

            for (i, bus) in data["bus"]
                if bus["bus_type"] != 4
                    index = parse(Int, i)
                    #println("$i, $(JuMP.value(opf_model[:va][index])), $(pm_sol["bus"][i]["va"])")
                    @test isapprox(JuMP.value(opf_model[:va][index]), pm_sol["bus"][i]["va"]; atol = 1e-8)
                end
            end

            for (i, gen) in data["gen"]
                if gen["gen_status"] != 0
                    index = parse(Int, i)
                    #println("$i, $(JuMP.value(opf_model[:pg][index])), $(pm_sol["gen"][i]["pg"])")
                    @test isapprox(JuMP.value(opf_model[:pg][index]), pm_sol["gen"][i]["pg"]; atol = 1e-8)
                end
            end
        else
            @test opf_status == LOCALLY_INFEASIBLE
            @test pm_result["status"] == LOCALLY_INFEASIBLE
        end
    end
end



function solve_file(file_name)
    include(file_name)
    optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
    return optimizer, data, termination_status(model), cost
end


@testset "test ac polar opf" begin
    optimizer, data, status, cost = solve_file("../../src/model/ac-opf.jl")
    pm_result = solve_ac_opf(data, optimizer)
    @test isapprox(cost, pm_result["objective"]; atol = 1e-6)
end


@testset "test dc polar opf" begin
    optimizer, data, status, cost = solve_file("../../src/model/dc-opf.jl")
    pm_result = solve_dc_opf(data, optimizer)
    @test isapprox(cost, pm_result["objective"]; atol = 1e-6)
end

