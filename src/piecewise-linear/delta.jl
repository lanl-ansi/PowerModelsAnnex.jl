export solve_opf_pwl_delta, build_opf_pwl_delta

""
function solve_opf_pwl_delta(file, model_type::Type, optimizer; kwargs...)
    return _PM.solve_model(file, model_type, optimizer, build_opf_pwl_delta; kwargs...)
end

"a variant of the OPF problem specification that uses a max variant of the pwl cost function implementation"
function build_opf_pwl_delta(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_dcline_power(pm)

    objective_min_fuel_cost_delta(pm)

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
function objective_min_fuel_cost_delta(pm::_PM.AbstractPowerModel; kwargs...)
    model = _PM.check_cost_models(pm)

    if model == 1
        return objective_min_fuel_cost_pwl_delta(pm; kwargs...)
    elseif model == 2
        return _PM.objective_min_fuel_cost_polynomial(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end


""
function objective_min_fuel_cost_pwl_delta(pm::_PM.AbstractPowerModel; kwargs...)
    objective_variable_pg_cost_delta(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_cost, i) for (i,gen) in nw_ref[:gen])
        for (n, nw_ref) in _PM.nws(pm))
    )
end


"adds pg_cost variables and constraints"
function objective_variable_pg_cost_delta(pm::_PM.AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in _PM.nws(pm)
        pg_cost = var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in ref(pm, n, :gen)
            points = _PM.calc_pwl_points(gen["ncost"], gen["cost"], gen["pmin"], gen["pmax"])

            cost_per_mw = Float64[0.0]
            for i in 2:length(points)
                x0 = points[i-1].mw
                y0 = points[i-1].cost
                x1 = points[i].mw
                y1 = points[i].cost

                m = (y1 - y0)/(x1 - x0)
                if !isnan(m)
                    push!(cost_per_mw, m)
                else
                    @assert isapprox(y0, y1)
                    push!(cost_per_mw, 0.0)
                end
            end

            pg_cost_mw = JuMP.@variable(pm.model,
                [i in 2:length(points)], base_name="$(n)_pg_cost_mw",
                lower_bound = 0.0,
                upper_bound = points[i].mw - points[i-1].mw
            )

            JuMP.@constraint(pm.model, points[1].mw + sum(pg_cost_mw[i] for i in 2:length(points)) == sum(var(pm, n, :pg, i)[c] for c in _PM.conductor_ids(pm, n)))
            pg_cost[i] = points[1].cost + sum(cost_per_mw[i]*pg_cost_mw[i] for i in 2:length(points))
        end

        report && _PM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)
    end
end
