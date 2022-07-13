export solve_opf_sad

""
function solve_opf_sad(file, model_constructor, optimizer; kwargs...)
    return _PM.solve_model(file, model_constructor, optimizer, build_opf_sad; kwargs...)
end

""
function build_opf_sad(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_dcline_power(pm, bounded = false)

    @variable(pm.model, theta_delta_bound >= 0.0, start = 0.523598776)

    @objective(pm.model, Min, theta_delta_bound)

    _PM.constraint_model_voltage(pm)

    for i in _PM.ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in _PM.ids(pm, :bus)
        _PM.constraint_power_balance(pm, i)
    end

    for (i, branch) in ref(pm, :branch)
        _PM.constraint_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)

        _PM.constraint_voltage_angle_difference(pm, i)
        theta_fr = var(pm, :va, branch["f_bus"])
        theta_to = var(pm, :va, branch["t_bus"])
        @constraint(pm.model, theta_fr - theta_to <=  theta_delta_bound)
        @constraint(pm.model, theta_fr - theta_to >= -theta_delta_bound)

        constraint_thermal_limit_from(pm, i; scale = 0.999)
        constraint_thermal_limit_to(pm, i; scale = 0.999)
    end

    for i in _PM.ids(pm, :dcline)
        _PM.constraint_dcline_power_losses(pm, i)
    end
end