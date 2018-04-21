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

    PMs.variable_branch_flow(pm)
    PMs.variable_dcline_flow(pm)

    variable_load_factor(pm)

    objective_max_loading(pm)
    #objective_max_loading_voltage_norm(pm)
    #objective_max_loading_gen_output(pm)

    PMs.constraint_voltage(pm)

    for i in ids(pm, :ref_buses)
        PMs.constraint_theta_ref(pm, i)
    end

    for (i,gen) in ref(pm, :gen)
        pg = var(pm, :pg, i)
        @constraint(pm.model, pg >= gen["pmin"])
    end

    for i in ids(pm, :bus)
        constraint_kcl_shunt_scaled(pm, i)
    end

    for i in ids(pm, :branch)
        PMs.constraint_ohms_yt_from(pm, i)
        PMs.constraint_ohms_yt_to(pm, i)

        PMs.constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i; scale = 0.999)
        constraint_thermal_limit_to(pm, i; scale = 0.999)
    end


    for i in ids(pm, :dcline)
        PMs.constraint_dcline(pm, i)
    end
end


"variable: load_factor >= 1.0"
function variable_load_factor(pm::GenericPowerModel)
    pm.var[:nw][pm.cnw][:load_factor] = @variable(pm.model,
        basename="load_factor",
        lowerbound=1.0,
        start = 1.0
    )
end

"objective: Max. load_factor"
function objective_max_loading(pm::GenericPowerModel)
    @objective(pm.model, Max, var(pm, :load_factor))
end

""
function objective_max_loading_voltage_norm(pm::GenericPowerModel)
    # Seems to create too much reactive power and makes even small models hard to converge
    load_factor = var(pm, :load_factor)

    scale = length(ids(pm, :bus))
    v = var(pm, :vm)

    @objective(pm.model, Max, 10*scale*load_factor - sum(((bus["vmin"] + bus["vmax"])/2 - v[i])^2 for (i,bus) in pm.ref[:bus]))
end

""
function objective_max_loading_gen_output(pm::GenericPowerModel)
    # Works but adds unnecessary runtime
    load_factor = var(pm, :load_factor)

    scale = length(ids(pm,:gen))
    pg = var(pm, :pg)
    qg = var(pm, :qg)

    @NLobjective(pm.model, Max, 100*scale*load_factor - sum( (pg[i]^2 - (2*qg[i])^2)^2 for (i,gen) in pm.ref[:gen] ))
end

""
function bounds_tighten_voltage{T <: PMs.AbstractACPForm}(pm::GenericPowerModel{T}; epsilon = 0.001)
    for (i,bus) in ref(pm, :bus)
        v = var(pm, :vm, i)
        setupperbound(v, bus["vmax"]*(1.0-epsilon))
        setlowerbound(v, bus["vmin"]*(1.0+epsilon))
    end
end

""
function upperbound_negative_active_generation(pm::GenericPowerModel)
    for (i,gen) in ref(pm, :gen)
        if gen["pmax"] <= 0
            pg = var(pm, :pg, i)
            setupperbound(pg, gen["pmax"])
        end
    end
end

""
function constraint_kcl_shunt_scaled{T <: PMs.AbstractACPForm}(pm::GenericPowerModel{T}, n::Int, i::Int)
    bus = pm.ref[:nw][n][:bus][i]
    bus_arcs = pm.ref[:nw][n][:bus_arcs][i]
    bus_gens = pm.ref[:nw][n][:bus_gens][i]
    bus_loads = pm.ref[:nw][n][:bus_loads][i]
    bus_shunts = pm.ref[:nw][n][:bus_shunts][i]

    load_factor = pm.var[:nw][n][:load_factor]
    v = pm.var[:nw][n][:vm]
    p = pm.var[:nw][n][:p]
    q = pm.var[:nw][n][:q]
    pg = pm.var[:nw][n][:pg]
    qg = pm.var[:nw][n][:qg]

    if length(bus_loads) > 0
        pd = sum([pm.ref[:nw][n][:load][i]["pd"] for i in bus_loads])
        qd = sum([pm.ref[:nw][n][:load][i]["qd"] for i in bus_loads])
    else
        pd = 0.0
        qd = 0.0
    end

    if length(bus_shunts) > 0 
        gs = sum([pm.ref[:nw][n][:shunt][i]["gs"] for i in bus_shunts])
        bs = sum([pm.ref[:nw][n][:shunt][i]["bs"] for i in bus_shunts])
    else
        gs = 0.0
        bs = 0.0
    end

    if pd > 0.0 && qd > 0.0
        @constraint(pm.model, sum(p[a] for a in bus_arcs) == sum(pg[g] for g in bus_gens) - pd*load_factor - gs*v[i]^2)
    else
        # super fallback impl
        @constraint(pm.model, sum(p[a] for a in bus_arcs) == sum(pg[g] for g in bus_gens) - pd - gs*v[i]^2)
    end

    @constraint(pm.model, sum(q[a] for a in bus_arcs) == sum(qg[g] for g in bus_gens) - qd + bs*v[i]^2)
end
constraint_kcl_shunt_scaled(pm::GenericPowerModel, i::Int) = constraint_kcl_shunt_scaled(pm, pm.cnw, i)


""
function get_api_solution(pm::GenericPowerModel, sol::Dict{String,Any})
    PMs.add_bus_voltage_setpoint(sol, pm)
    PMs.add_generator_power_setpoint(sol, pm)
    PMs.add_branch_flow_setpoint(sol, pm)

    # extension
    add_load_demand_setpoint(sol, pm)
end

""
function add_load_demand_setpoint(sol, pm::GenericPowerModel)
    mva_base = pm.data["baseMVA"]
    PMs.add_setpoint(sol, pm, "load", "pd", :load_factor; default_value = (item) -> item["pd"], scale = (x,item) -> item["pd"] > 0 && item["qd"] > 0 ? x*item["pd"] : item["pd"], extract_var = (var,idx,item) -> var)
    PMs.add_setpoint(sol, pm, "load", "qd", :load_factor; default_value = (item) -> item["qd"], scale = (x,item) -> item["qd"], extract_var = (var,idx,item) -> var)
end
