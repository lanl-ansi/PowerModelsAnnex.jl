export run_api_opf

""
function run_api_opf(file, model_constructor, optimizer; kwargs...)
    return PMs.run_model(file, model_constructor, optimizer, post_api_opf; solution_builder = solution_api, kwargs...)
end

""
function post_api_opf(pm::PMs.GenericPowerModel)
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
function variable_load_factor(pm::PMs.GenericPowerModel)
    var(pm)[:load_factor] = @variable(pm.model,
        base_name="load_factor",
        lower_bound=1.0,
        start = 1.0
    )
end

"objective: Max. load_factor"
function objective_max_loading(pm::PMs.GenericPowerModel)
    @objective(pm.model, Max, var(pm, :load_factor))
end

""
function objective_max_loading_voltage_norm(pm::PMs.GenericPowerModel)
    # Seems to create too much reactive power and makes even small models hard to converge
    load_factor = var(pm, :load_factor)

    scale = length(PMs.ids(pm, :bus))
    vm = var(pm, :vm)

    @objective(pm.model, Max, 10*scale*load_factor - sum(((bus["vmin"] + bus["vmax"])/2 - vm[i])^2 for (i,bus) in ref(pm, :bus)))
end

""
function objective_max_loading_gen_output(pm::PMs.GenericPowerModel)
    # Works but adds unnecessary runtime
    load_factor = var(pm, :load_factor)

    scale = length(PMs.ids(pm, :gen))
    pg = var(pm, :pg)
    qg = var(pm, :qg)

    @NLobjective(pm.model, Max, 100*scale*load_factor - sum( (pg[i]^2 - (2*qg[i])^2)^2 for (i,gen) in ref(pm, :gen)))
end

""
function bounds_tighten_voltage(pm::PMs.GenericPowerModel{T}; epsilon = 0.001) where T <: PMs.AbstractACPForm
    for (i,bus) in ref(pm, :bus)
        v = var(pm, :vm, i)
        JuMP.set_upper_bound(v, bus["vmax"]*(1.0-epsilon))
        JuMP.set_lower_bound(v, bus["vmin"]*(1.0+epsilon))
    end
end

""
function upperbound_negative_active_generation(pm::PMs.GenericPowerModel)
    for (i,gen) in ref(pm, :gen)
        if gen["pmax"] <= 0
            pg = var(pm, :pg, i)
            JuMP.set_upper_bound(pg, gen["pmax"])
        end
    end
end

""
function constraint_power_balance_shunt_scaled(pm::PMs.GenericPowerModel{T}, n::Int, c::Int, i::Int) where T <: PMs.AbstractACPForm
    bus = ref(pm, n, :bus, i)
    bus_arcs = ref(pm, n, :bus_arcs, i)
    bus_gens = ref(pm, n, :bus_gens, i)
    bus_loads = ref(pm, n, :bus_loads, i)
    bus_shunts = ref(pm, n, :bus_shunts, i)

    load_factor = var(pm, n, c, :load_factor)
    vm = var(pm, n, c, :vm, i)
    p = var(pm, n, c, :p)
    q = var(pm, n, c, :q)
    pg = var(pm, n, c, :pg)
    qg = var(pm, n, c, :qg)

    if length(bus_loads) > 0
        pd = sum([ref(pm, n, :load, i, "pd", c) for i in bus_loads])
        qd = sum([ref(pm, n, :load, i, "qd", c) for i in bus_loads])
    else
        pd = 0.0
        qd = 0.0
    end

    if length(bus_shunts) > 0
        gs = sum([ref(pm, n, :shunt, i, "gs", c) for i in bus_shunts])
        bs = sum([ref(pm, n, :shunt, i, "bs", c) for i in bus_shunts])
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
constraint_power_balance_shunt_scaled(pm::PMs.GenericPowerModel, i::Int) = constraint_power_balance_shunt_scaled(pm, pm.cnw, pm.ccnd, i)


""
function solution_api(pm::PMs.GenericPowerModel, sol::Dict{String,Any})
    PMs.add_setpoint_bus_voltage!(sol, pm)
    PMs.add_setpoint_generator_power!(sol, pm)
    PMs.add_setpoint_branch_flow!(sol, pm)

    # extension
    add_setpoint_load_demand!(sol, pm)
end

Base.getindex(x::JuMP.VariableRef, i::Int64) = x

""
function add_setpoint_load_demand!(sol, pm::PMs.GenericPowerModel)
    mva_base = pm.data["baseMVA"]
    PMs.add_setpoint!(sol, pm, "load", "pd", :load_factor; default_value = (item) -> item["pd"], scale = (x,item,i) -> item["pd"][i] > 0 && item["qd"][i] > 0 ? x*item["pd"][i] : item["pd"][i], var_key = (idx,item) -> 1)
    PMs.add_setpoint!(sol, pm, "load", "qd", :load_factor; default_value = (item) -> item["qd"], scale = (x,item,i) -> item["qd"][i], var_key = (idx,item) -> 1)
end