### Branch - Thermal Limit Constraints ###

""
function constraint_thermal_limit_from(pm::GenericPowerModel, n::Int, c::Int, i::Int; scale = 1.0)
    branch = ref(pm, n, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    PMs.constraint_thermal_limit_from(pm, n, c, f_idx, branch["rate_a"][c]*scale)
end
constraint_thermal_limit_from(pm::GenericPowerModel, i::Int; scale = 1.0) = constraint_thermal_limit_from(pm, pm.cnw, pm.ccnd, i; scale=scale)

""
function constraint_thermal_limit_to(pm::GenericPowerModel, n::Int, c::Int, i::Int; scale = 1.0)
    branch = ref(pm, n, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    PMs.constraint_thermal_limit_to(pm, n, c, t_idx, branch["rate_a"][c]*scale)
end
constraint_thermal_limit_to(pm::GenericPowerModel, i::Int; scale = 1.0) = constraint_thermal_limit_to(pm, pm.cnw, pm.ccnd, i; scale=scale)


