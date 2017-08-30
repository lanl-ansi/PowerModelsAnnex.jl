export run_api_opf

""
function run_api_opf(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_api_opf; solution_builder = get_api_solution, kwargs...)
end

""
function post_api_opf(pm::GenericPowerModel)
    PMs.variable_voltage(pm)
    bounds_tighten_voltage(pm)

    PMs.variable_generation(pm, bounded = false)
    upperbound_negative_active_generation(pm)

    PMs.variable_line_flow(pm)
    PMs.variable_dcline_flow(pm)

    variable_load_factor(pm)

    objective_max_loading(pm)
    #objective_max_loading_voltage_norm(pm)
    #objective_max_loading_gen_output(pm)

    PMs.constraint_voltage(pm)

    for (i,bus) in pm.ref[:ref_buses]
        PMs.constraint_theta_ref(pm, bus)
    end

    for (i,gen) in pm.ref[:gen]
        pg = pm.var[:pg][i]
        @constraint(pm.model, pg >= gen["pmin"])
    end

    for (i,bus) in pm.ref[:bus]
        constraint_kcl_shunt_scaled(pm, bus)
    end

    for (i,branch) in pm.ref[:branch]
        PMs.constraint_ohms_yt_from(pm, branch)
        PMs.constraint_ohms_yt_to(pm, branch)

        PMs.constraint_voltage_angle_difference(pm, branch)

        constraint_thermal_limit_from(pm, branch; scale = 0.999)
        constraint_thermal_limit_to(pm, branch; scale = 0.999)
    end


    for (i,dcline) in pm.ref[:dcline]
        PMs.constraint_dcline(pm, dcline)
    end
end


"variable: load_factor >= 1.0"
function variable_load_factor(pm::GenericPowerModel)
    pm.var[:load_factor] = @variable(pm.model,
        basename="load_factor",
        lowerbound=1.0,
        start = 1.0
    )
end

"objective: Max. load_factor"
function objective_max_loading(pm::GenericPowerModel)
    @objective(pm.model, Max, pm.var[:load_factor])
end

""
function objective_max_loading_voltage_norm(pm::GenericPowerModel)
    # Seems to create too much reactive power and makes even small models hard to converge
    load_factor = pm.var[:load_factor]

    scale = length(pm.ref[:bus])
    v = pm.var[:vm]

    @objective(pm.model, Max, 10*scale*load_factor - sum(((bus["vmin"] + bus["vmax"])/2 - v[i])^2 for (i,bus) in pm.ref[:bus]))
end

""
function objective_max_loading_gen_output(pm::GenericPowerModel)
    # Works but adds unnecessary runtime
    load_factor = pm.var[:load_factor]

    scale = length(pm.ref[:gen])
    pg = pm.var[:pg]
    qg = pm.var[:qg]

    @NLobjective(pm.model, Max, 100*scale*load_factor - sum( (pg[i]^2 - (2*qg[i])^2)^2 for (i,gen) in pm.ref[:gen] ))
end

""
function bounds_tighten_voltage{T <: PMs.AbstractACPForm}(pm::GenericPowerModel{T}; epsilon = 0.001)
    for (i,bus) in pm.ref[:bus]
        v = pm.var[:vm][i]
        setupperbound(v, bus["vmax"]*(1.0-epsilon))
        setlowerbound(v, bus["vmin"]*(1.0+epsilon))
    end
end

""
function upperbound_negative_active_generation(pm::GenericPowerModel)
    for (i,gen) in pm.ref[:gen]
        if gen["pmax"] <= 0
            pg = pm.var[:pg][i]
            setupperbound(pg, gen["pmax"])
        end
    end
end

""
function constraint_kcl_shunt_scaled{T <: PMs.AbstractACPForm}(pm::GenericPowerModel{T}, bus)
    i = bus["index"]
    bus_arcs = pm.ref[:bus_arcs][i]
    bus_gens = pm.ref[:bus_gens][i]

    load_factor = pm.var[:load_factor]
    v = pm.var[:vm]
    p = pm.var[:p]
    q = pm.var[:q]
    pg = pm.var[:pg]
    qg = pm.var[:qg]

    if bus["pd"] > 0 && bus["qd"] > 0
        @constraint(pm.model, sum(p[a] for a in bus_arcs) == sum(pg[g] for g in bus_gens) - bus["pd"]*load_factor - bus["gs"]*v[i]^2)
    else
        # super fallback impl
        @constraint(pm.model, sum(p[a] for a in bus_arcs) == sum(pg[g] for g in bus_gens) - bus["pd"] - bus["gs"]*v[i]^2)
    end

    @constraint(pm.model, sum(q[a] for a in bus_arcs) == sum(qg[g] for g in bus_gens) - bus["qd"] + bus["bs"]*v[i]^2)
end

""
function get_api_solution(pm::GenericPowerModel)
    # super fallback
    sol = PMs.init_solution(pm)
    PMs.add_bus_voltage_setpoint(sol, pm)
    PMs.add_generator_power_setpoint(sol, pm)
    PMs.add_branch_flow_setpoint(sol, pm)

    # extension
    add_bus_demand_setpoint(sol, pm)

    return sol
end

""
function add_bus_demand_setpoint(sol, pm::GenericPowerModel)
    mva_base = pm.data["baseMVA"]
    PMs.add_setpoint(sol, pm, "bus", "pd", :load_factor; default_value = (item) -> item["pd"], scale = (x,item) -> item["pd"] > 0 && item["qd"] > 0 ? x*item["pd"] : item["pd"], extract_var = (var,idx,item) -> var)
    PMs.add_setpoint(sol, pm, "bus", "qd", :load_factor; default_value = (item) -> item["qd"], scale = (x,item) -> item["qd"], extract_var = (var,idx,item) -> var)
end
