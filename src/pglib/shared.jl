### Branch - Thermal Limit Constraints ###

""
function constraint_thermal_limit_from(pm::PMs.AbstractPowerModel, n::Int, i::Int; scale = 1.0)
    branch = ref(pm, n, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if haskey(branch, "rate_a")
        PMs.constraint_thermal_limit_from(pm, n, f_idx, branch["rate_a"]*scale)
    end
end
constraint_thermal_limit_from(pm::PMs.AbstractPowerModel, i::Int; scale = 1.0) = constraint_thermal_limit_from(pm, pm.cnw, i; scale=scale)

""
function constraint_thermal_limit_to(pm::PMs.AbstractPowerModel, n::Int, i::Int; scale = 1.0)
    branch = ref(pm, n, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    if haskey(branch, "rate_a")
        PMs.constraint_thermal_limit_to(pm, n, t_idx, branch["rate_a"]*scale)
    end
end
constraint_thermal_limit_to(pm::PMs.AbstractPowerModel, i::Int; scale = 1.0) = constraint_thermal_limit_to(pm, pm.cnw, i; scale=scale)


