export run_opf_pwl_psi, build_opf_pwl_psi

""
function run_opf_pwl_psi(file, model_type::Type, optimizer; kwargs...)
    return _PM.run_model(file, model_type, optimizer, build_opf_pwl_psi; kwargs...)
end

"a variant of the OPF problem specification that uses a max variant of the pwl cost function implementation"
function build_opf_pwl_psi(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_psi(pm)

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

""
function objective_min_fuel_and_flow_cost_psi(pm::_PM.AbstractPowerModel; kwargs...)
    model = _PM.check_cost_models(pm)

    if model == 1
        return objective_min_fuel_and_flow_cost_pwl_psi(pm; kwargs...)
    elseif model == 2
        return objective_min_fuel_and_flow_cost_polynomial(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end

""
function objective_min_fuel_and_flow_cost_pwl_psi(pm::_PM.AbstractPowerModel; kwargs...)
    objective_variable_pg_cost_psi(pm; kwargs...)
    _PM.objective_variable_dc_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_cost, i) for (i,gen) in nw_ref[:gen]) +
            sum( var(pm, n, :p_dc_cost, i) for (i,dcline) in nw_ref[:dcline])
        for (n, nw_ref) in _PM.nws(pm))
    )
end

"adds pg_cost variables and constraints"
function objective_variable_pg_cost_psi(pm::_PM.AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in _PM.nws(pm)
        gen_lines = Dict{Int64,Vector}()
        pg_cost_start = Dict{Int64,Float64}()
        pg_cost_min = Dict{Int64,Float64}()
        pg_cost_max = Dict{Int64,Float64}()

        for (i, gen) in nw_ref[:gen]
            @assert gen["model"] == 1
            ncost, points = get_active_cost_points(gen)

            mws = Float64[]
            costs = Float64[]
            for i in 1:ncost
                push!(mws, points[2*i-1])
                push!(costs, points[2*i])
            end

            gen_lines[i] = calc_comp_lines_from_points(points)
            pg_value = sum(JuMP.start_value(var(pm, n, :pg, i)[c]) for c in _PM.conductor_ids(pm, n))
            pg_cost_value = -Inf
            for line in gen_lines[i]
                pg_cost_value = max(pg_cost_value, line.slope*pg_value + line.intercept)
            end
            pg_cost_start[i] = pg_cost_value
            pg_cost_min[i] = costs[1]
            pg_cost_max[i] = costs[end]
        end


        pg_cost = var(pm, n)[:pg_cost] = JuMP.@variable(pm.model,
            [i in ids(pm, n, :gen)], base_name="$(n)_pg_cost",
            lower_bound = pg_cost_min[i],
            upper_bound = pg_cost_max[i]
        )
        report && _IM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)

        # gen pwl cost
        for (i, gen) in nw_ref[:gen]
            for line in gen_lines[i]
                JuMP.@constraint(pm.model, pg_cost[i] >= line.slope*sum(var(pm, n, :pg, i)[c] for c in _PM.conductor_ids(pm, n)) + line.intercept)
            end
        end
    end
end

