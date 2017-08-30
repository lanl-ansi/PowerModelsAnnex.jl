export run_sad_opf

""
function run_sad_opf(file, model_constructor, solver; kwargs...)
    return run_generic_model(file, model_constructor, solver, post_sad_opf; kwargs...)
end

""
function post_sad_opf{T}(pm::GenericPowerModel{T})
    PMs.variable_voltage(pm)
    PMs.variable_generation(pm)
    PMs.variable_line_flow(pm)
    PMs.variable_dcline_flow(pm, bounded = false)

    @variable(pm.model, theta_delta_bound >= 0.0, start = 0.523598776)

    @objective(pm.model, Min, theta_delta_bound)

    PMs.constraint_voltage(pm)

    for (i,bus) in pm.ref[:ref_buses]
        PMs.constraint_theta_ref(pm, bus)
    end

    for (i,bus) in pm.ref[:bus]
        PMs.constraint_kcl_shunt(pm, bus)
    end

    for (i,branch) in pm.ref[:branch]
        PMs.constraint_ohms_yt_from(pm, branch)
        PMs.constraint_ohms_yt_to(pm, branch)

        PMs.constraint_voltage_angle_difference(pm, branch)
        theta_fr = pm.var[:va][branch["f_bus"]]
        theta_to = pm.var[:va][branch["t_bus"]]

        @constraint(pm.model, theta_fr - theta_to <=  theta_delta_bound)
        @constraint(pm.model, theta_fr - theta_to >= -theta_delta_bound)

        constraint_thermal_limit_from(pm, branch; scale = 0.999)
        constraint_thermal_limit_to(pm, branch; scale = 0.999)
    end

    for (i,dcline) in pm.ref[:dcline]
        PMs.constraint_dcline(pm, dcline)
    end
end
