export solve_opf_cop

"a closet operating point (cop) opf model variant that minimizes distance from the given operating point"
function solve_opf_cop(file, model_type::Type, optimizer; kwargs...)
    return _PM.solve_model(file, model_type, optimizer, build_opf_cop; kwargs...)
end

""
function build_opf_cop(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_dcline_power(pm)

    objective_min_pg_vm_distance(pm)

    _PM.constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        _PM.constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        _PM.constraint_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)

        _PM.constraint_voltage_angle_difference(pm, i)

        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        _PM.constraint_dcline_power_losses(pm, i)
    end
end


function objective_min_pg_vm_distance(pm::_PM.AbstractPowerModel)
    nws = _PM.nw_ids(pm)
    @assert all(!_PM.ismulticonductor(pm, n) for n in nws)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum((gen["pg"] - var(pm, n, :pg, i))^2 for (i,gen) in ref(pm, n, :gen)) +
            sum((bus["vm"] - var(pm, n, :vm, i))^2 for (i,bus) in ref(pm, n, :bus)) 
        for n in nws)
    )
end

function objective_min_pg_vm_distance(pm::_PM.AbstractActivePowerModel)
    nws = _PM.nw_ids(pm)
    @assert all(!_PM.ismulticonductor(pm, n) for n in nws)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum((gen["pg"] - var(pm, n, :pg, i))^2 for (i,gen) in ref(pm, n, :gen))
        for n in nws)
    )
end