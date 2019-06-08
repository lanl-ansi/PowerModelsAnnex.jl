export run_sad_opf

""
function run_sad_opf(file, model_constructor, optimizer; kwargs...)
    return PMs.run_model(file, model_constructor, optimizer, post_sad_opf; kwargs...)
end

""
function post_sad_opf(pm::PMs.GenericPowerModel)
    PMs.variable_voltage(pm)
    PMs.variable_generation(pm)
    PMs.variable_branch_flow(pm)
    PMs.variable_dcline_flow(pm, bounded = false)

    @variable(pm.model, theta_delta_bound >= 0.0, start = 0.523598776)

    @objective(pm.model, Min, theta_delta_bound)

    PMs.constraint_model_voltage(pm)

    for i in PMs.ids(pm, :ref_buses)
        PMs.constraint_theta_ref(pm, i)
    end

    for i in PMs.ids(pm, :bus)
        PMs.constraint_power_balance_shunt(pm, i)
    end

    for (i, branch) in ref(pm, :branch)
        PMs.constraint_ohms_yt_from(pm, i)
        PMs.constraint_ohms_yt_to(pm, i)

        PMs.constraint_voltage_angle_difference(pm, i)
        theta_fr = var(pm, :va, branch["f_bus"])
        theta_to = var(pm, :va, branch["t_bus"])
        @constraint(pm.model, theta_fr - theta_to <=  theta_delta_bound)
        @constraint(pm.model, theta_fr - theta_to >= -theta_delta_bound)

        constraint_thermal_limit_from(pm, i; scale = 0.999)
        constraint_thermal_limit_to(pm, i; scale = 0.999)
    end

    for i in PMs.ids(pm, :dcline)
        PMs.constraint_dcline(pm, i)
    end
end