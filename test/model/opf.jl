
function run_ac_opf_model(data, solver)
    model = post_ac_opf(data, JuMP.Model(solver))
    JuMP.optimize!(model)
    status = PMs.parse_status(JuMP.termination_status(model), JuMP.primal_status(model), JuMP.dual_status(model))
    return status, model
end

@testset "test ac polar opf" begin
    @testset "case $(name)" for (name, case_file) in case_files
        data = PMs.parse_file(case_file)
        opf_status, opf_model = run_ac_opf_model(data, ipopt_solver)
        pm_result = PMs.run_ac_opf(data, ipopt_solver)
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

function run_soc_opf_model(data, solver)
    model =  post_soc_opf(data, JuMP.Model(solver))
    JuMP.optimize!(model)
    status = PMs.parse_status(JuMP.termination_status(model), JuMP.primal_status(model), JuMP.dual_status(model))
    return status, model
end

@testset "test soc w opf" begin
    @testset "case $(name)" for (name, case_file) in case_files
        data = PMs.parse_file(case_file)
        opf_status, opf_model = run_soc_opf_model(data, ipopt_solver)
        pm_result = PMs.run_opf(data, PMs.SOCWRPowerModel, ipopt_solver)
        pm_sol = pm_result["solution"]

        @test isapprox(JuMP.objective_value(opf_model), pm_result["objective"]; atol = 1e-5)

        base_mva = data["baseMVA"]

        for (i, bus) in data["bus"]
            if bus["bus_type"] != 4
                index = parse(Int, i)
                #println("$i, $(JuMP.value(opf_model[:va][index])), $(pm_sol["bus"][i]["va"])")
                #@test isapprox(JuMP.value(opf_model[:va][index]), pm_sol["bus"][i]["va"]; atol = 1e-8)

                #println("$i, $(JuMP.value(opf_model[:w][index])), $(pm_sol["bus"][i]["vm"]^2)")
                @test isapprox(JuMP.value(opf_model[:w][index]), pm_sol["bus"][i]["vm"]^2; atol = 1e-6)
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


function run_qc_opf_model(data, solver)
    model =  post_qc_opf(data, JuMP.Model(solver))
    JuMP.optimize!(model)
    status = PMs.parse_status(JuMP.termination_status(model), JuMP.primal_status(model), JuMP.dual_status(model))
    return status, model
end

@testset "test qc w+l opf" begin
    @testset "case $(name)" for (name, case_file) in case_files
        data = PMs.parse_file(case_file)
        opf_status, opf_model = run_qc_opf_model(data, ipopt_solver)
        pm_result = PMs.run_opf(data, PMs.QCWRTriPowerModel, ipopt_solver)
        pm_sol = pm_result["solution"]

        @test isapprox(JuMP.objective_value(opf_model), pm_result["objective"]; atol = 1e-5)

        base_mva = data["baseMVA"]

        for (i, bus) in data["bus"]
            if bus["bus_type"] != 4
                index = parse(Int, i)
                #println("$i, $(JuMP.value(opf_model[:va][index])), $(pm_sol["bus"][i]["va"])")
                #@test isapprox(JuMP.value(opf_model[:va][index]), pm_sol["bus"][i]["va"]; atol = 1e-8)

                #println("$i, $(JuMP.value(opf_model[:vm][index])), $(pm_sol["bus"][i]["vm"])")
                @test isapprox(JuMP.value(opf_model[:vm][index]), pm_sol["bus"][i]["vm"]; atol = 1e-6)
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


function run_dc_opf_model(data, solver)
    model = post_dc_opf(data, JuMP.Model(solver))
    JuMP.optimize!(model)
    status = PMs.parse_status(JuMP.termination_status(model), JuMP.primal_status(model), JuMP.dual_status(model))
    return status, model
end

@testset "test dc polar opf" begin
    @testset "case $(name)" for (name, case_file) in case_files
        data = PMs.parse_file(case_file)
        opf_status, opf_model = run_dc_opf_model(data, ipopt_solver)
        pm_result = PMs.run_dc_opf(data, ipopt_solver)
        pm_sol = pm_result["solution"]

        #println(opf_status)
        #println(pm_result["status"])

        @test isapprox(JuMP.objective_value(opf_model), pm_result["objective"]; atol = 1e-5)

        # needed becouse some test networks are not DC feasible
        if pm_result["status"] == :LocalOptimal
            @test opf_status == :LocalOptimal

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
                    @test isapprox(JuMP.value(opf_model[:pg][index]), pm_sol["gen"][i]["pg"])
                end
            end
        else
            @test opf_status == :Infeasible
            @test pm_result["status"] == :LocalInfeasible
        end
    end
end



function run_file(file_name)
    include(file_name)
    return nlp_solver, data, status, cost
end


@testset "test ac polar opf" begin
    solver, data, status, cost = run_file("../../src/model/ac-opf.jl")
    pm_result = PMs.run_ac_opf(data, solver)
    @test isapprox(cost, pm_result["objective"]; atol = 1e-6)
end


@testset "test dc polar opf" begin
    solver, data, status, cost = run_file("../../src/model/dc-opf.jl")
    pm_result = PMs.run_dc_opf(data, solver)
    @test isapprox(cost, pm_result["objective"]; atol = 1e-6)
end

