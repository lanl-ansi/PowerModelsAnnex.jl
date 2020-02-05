export run_opf_api

""
function run_opf_api(file, model_constructor, optimizer; kwargs...)
    return PMs.run_model(file, model_constructor, optimizer, build_opf_api; kwargs...)
end

""
function build_opf_api(pm::PMs.AbstractPowerModel)
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

    PMs.constraint_model_voltage(pm)

    for i in PMs.ids(pm, :ref_buses)
        PMs.constraint_theta_ref(pm, i)
    end

    for (i,gen) in ref(pm, :gen)
        pg = var(pm, :pg, i)
        @constraint(pm.model, pg >= gen["pmin"])
    end

    for i in PMs.ids(pm, :bus)
        constraint_power_balance_shunt_scaled(pm, i)
    end

    for i in PMs.ids(pm, :branch)
        PMs.constraint_ohms_yt_from(pm, i)
        PMs.constraint_ohms_yt_to(pm, i)

        PMs.constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i; scale = 0.999)
        constraint_thermal_limit_to(pm, i; scale = 0.999)
    end


    for i in PMs.ids(pm, :dcline)
        PMs.constraint_dcline(pm, i)
    end
end


"variable: load_factor >= 1.0"
function variable_load_factor(pm::PMs.AbstractPowerModel, report::Bool=true)
    load_factor = var(pm)[:load_factor] = @variable(pm.model,
        base_name="load_factor",
        lower_bound=1.0,
        start = 1.0
    )
    sol(pm)[:load_factor] = load_factor
    for (i,load) in ref(pm, :load)
        #item["pd"][i] > 0 && item["qd"][i] > 0 ? x*item["pd"][i] : item["pd"][i]
        if load["pd"] > 0 && load["qd"] > 0
            sol(pm, :load, i)[:pd] = load["pd"]*load_factor
        else
            sol(pm, :load, i)[:pd] = load["pd"]
        end
        sol(pm, :load, i)[:qd] = load["qd"]*load_factor
    end
    #PMs.add_setpoint!(sol, pm, "load", "pd", :load_factor; default_value = (item) -> item["pd"], scale = (x,item,i) -> item["pd"][i] > 0 && item["qd"][i] > 0 ? x*item["pd"][i] : item["pd"][i], var_key = (idx,item) -> 1)
    #PMs.add_setpoint!(sol, pm, "load", "qd", :load_factor; default_value = (item) -> item["qd"], scale = (x,item,i) -> item["qd"][i], var_key = (idx,item) -> 1)
end

"objective: Max. load_factor"
function objective_max_loading(pm::PMs.AbstractPowerModel)
    @objective(pm.model, Max, var(pm, :load_factor))
end

""
function objective_max_loading_voltage_norm(pm::PMs.AbstractPowerModel)
    # Seems to create too much reactive power and makes even small models hard to converge
    load_factor = var(pm, :load_factor)

    scale = length(PMs.ids(pm, :bus))
    vm = var(pm, :vm)

    @objective(pm.model, Max, 10*scale*load_factor - sum(((bus["vmin"] + bus["vmax"])/2 - vm[i])^2 for (i,bus) in ref(pm, :bus)))
end

""
function objective_max_loading_gen_output(pm::PMs.AbstractPowerModel)
    # Works but adds unnecessary runtime
    load_factor = var(pm, :load_factor)

    scale = length(PMs.ids(pm, :gen))
    pg = var(pm, :pg)
    qg = var(pm, :qg)

    @NLobjective(pm.model, Max, 100*scale*load_factor - sum( (pg[i]^2 - (2*qg[i])^2)^2 for (i,gen) in ref(pm, :gen)))
end

""
function bounds_tighten_voltage(pm::PMs.AbstractACPModel; epsilon = 0.001)
    for (i,bus) in ref(pm, :bus)
        v = var(pm, :vm, i)
        JuMP.set_upper_bound(v, bus["vmax"]*(1.0-epsilon))
        JuMP.set_lower_bound(v, bus["vmin"]*(1.0+epsilon))
    end
end

""
function upperbound_negative_active_generation(pm::PMs.AbstractPowerModel)
    for (i,gen) in ref(pm, :gen)
        if gen["pmax"] <= 0
            pg = var(pm, :pg, i)
            JuMP.set_upper_bound(pg, gen["pmax"])
        end
    end
end

""
function constraint_power_balance_shunt_scaled(pm::PMs.AbstractACPModel, n::Int, i::Int)
    bus = ref(pm, n, :bus, i)
    bus_arcs = ref(pm, n, :bus_arcs, i)
    bus_gens = ref(pm, n, :bus_gens, i)
    bus_loads = ref(pm, n, :bus_loads, i)
    bus_shunts = ref(pm, n, :bus_shunts, i)

    load_factor = var(pm, n, :load_factor)
    vm = var(pm, n, :vm, i)
    p = var(pm, n, :p)
    q = var(pm, n, :q)
    pg = var(pm, n, :pg)
    qg = var(pm, n, :qg)

    if length(bus_loads) > 0
        pd = sum([ref(pm, n, :load, i, "pd") for i in bus_loads])
        qd = sum([ref(pm, n, :load, i, "qd") for i in bus_loads])
    else
        pd = 0.0
        qd = 0.0
    end

    if length(bus_shunts) > 0
        gs = sum([ref(pm, n, :shunt, i, "gs") for i in bus_shunts])
        bs = sum([ref(pm, n, :shunt, i, "bs") for i in bus_shunts])
    else
        gs = 0.0
        bs = 0.0
    end

    if pd > 0.0 && qd > 0.0
        @constraint(pm.model, sum(p[a] for a in bus_arcs) == sum(pg[g] for g in bus_gens) - pd*load_factor - gs*vm^2)
    else
        # super fallback impl
        @constraint(pm.model, sum(p[a] for a in bus_arcs) == sum(pg[g] for g in bus_gens) - pd - gs*vm^2)
    end

    @constraint(pm.model, sum(q[a] for a in bus_arcs) == sum(qg[g] for g in bus_gens) - qd + bs*vm^2)
end
constraint_power_balance_shunt_scaled(pm::PMs.AbstractPowerModel, i::Int) = constraint_power_balance_shunt_scaled(pm, pm.cnw, i)

