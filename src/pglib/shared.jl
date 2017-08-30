""
function constraint_thermal_limit_from(pm::GenericPowerModel, branch; scale = 1.0)
    i = branch["index"]
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    PMs.constraint_thermal_limit_from(pm, f_idx, branch["rate_a"]*scale)
end

""
function constraint_thermal_limit_to(pm::GenericPowerModel, branch; scale = 1.0)
    i = branch["index"]
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    PMs.constraint_thermal_limit_to(pm, t_idx, branch["rate_a"]*scale)
end
