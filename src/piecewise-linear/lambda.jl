export solve_opf_pwl_lambda, build_opf_pwl_lambda

""
function solve_opf_pwl_lambda(file, model_type::Type, optimizer; kwargs...)
    return _PM.solve_model(file, model_type, optimizer, build_opf_pwl_lambda; kwargs...)
end

"a variant of the OPF problem specification that uses a max variant of the pwl cost function implementation"
function build_opf_pwl_lambda(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_dcline_power(pm)

    objective_min_fuel_cost_lambda(pm)

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
function objective_min_fuel_cost_lambda(pm::_PM.AbstractPowerModel; kwargs...)
    expression_pg_cost_lambda(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n, :pg_cost, i) for (i,gen) in nw_ref[:gen])
        for (n, nw_ref) in _PM.nws(pm))
    )
end

""
function expression_pg_cost_lambda(pm::_PM.AbstractPowerModel; report::Bool=true)
    for (n, nw_ref) in _PM.nws(pm)
        pg_cost = var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in ref(pm, n, :gen)
            pg_terms = [var(pm, n, :pg, i)]

            if gen["model"] == 1
                if isa(pg_terms, Array{JuMP.VariableRef})
                    pmin = sum(JuMP.lower_bound.(pg_terms))
                    pmax = sum(JuMP.upper_bound.(pg_terms))
                else
                    pmin = gen["pmin"]
                    pmax = gen["pmax"]
                end

                points = _PM.calc_pwl_points(gen["ncost"], gen["cost"], pmin, pmax)
                pg_cost[i] = _pwl_cost_expression_lambda(pm, pg_terms, points, nw=n, id=i, var_name="pg")

            elseif gen["model"] == 2
                cost_rev = reverse(gen["cost"])

                pg_cost[i] = _PM._polynomial_cost_expression(pm, pg_terms, cost_rev, nw=n, id=i, var_name="pg")
            else
                Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model) on generator $(i)")
            end
        end

        report && _PM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)
    end
end

""
function _pwl_cost_expression_lambda(pm::_PM.AbstractPowerModel, x_list::Array{JuMP.VariableRef}, points; nw=0, id=1, var_name="x")
    pg_cost_lambda = JuMP.@variable(pm.model,
        [i in 1:length(points)], base_name="$(nw)_pg_$(id)_cost_lambda_$(i)",
        lower_bound = 0.0,
        upper_bound = 1.0
    )
    JuMP.@constraint(pm.model, sum(pg_cost_lambda) == 1.0)

    pg_expr = sum(pt.mw*pg_cost_lambda[i] for (i,pt) in enumerate(points))
    pg_cost_expr = sum(pt.cost*pg_cost_lambda[i] for (i,pt) in enumerate(points))

    JuMP.@constraint(pm.model, pg_expr == var(pm, nw, :pg, id))

    return pg_cost_expr
end

