export solve_opf_pwl_psi, build_opf_pwl_psi

""
function solve_opf_pwl_psi(file, model_type::Type, optimizer; kwargs...)
    return _PM.solve_model(file, model_type, optimizer, build_opf_pwl_psi; kwargs...)
end

"a variant of the OPF problem specification that uses a max variant of the pwl cost function implementation"
function build_opf_pwl_psi(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_dcline_power(pm)

    objective_min_fuel_cost_psi(pm)

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
function objective_min_fuel_cost_psi(pm::_PM.AbstractPowerModel; kwargs...)
    nl_gen = _PM.check_nl_gen_cost_models(pm)

    nl = nl_gen || typeof(pm) <: _PM.AbstractIVRModel

    expression_pg_cost_psi(pm; kwargs...)

    if !nl
        return JuMP.@objective(pm.model, Min,
            sum(
                sum( var(pm, n, :pg_cost, i) for (i,gen) in nw_ref[:gen])
            for (n, nw_ref) in _PM.nws(pm))
        )
    else
        pg_cost = Dict()
        for (n, nw_ref) in nws(pm)
            for (i,gen) in nw_ref[:gen]
                pg_cost[(n,i)] = var(pm, n, :pg_cost, i)
            end
        end

        return JuMP.@NLobjective(pm.model, Min,
            sum(
                sum( pg_cost[n,i] for (i,gen) in nw_ref[:gen])
            for (n, nw_ref) in _PM.nws(pm))
        )
    end
end

""
function expression_pg_cost_psi(pm::_PM.AbstractPowerModel; report::Bool=true)
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
                pg_cost[i] = _pwl_cost_expression_psi(pm, pg_terms, points, nw=n, id=i, var_name="pg")

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
function _pwl_cost_expression_psi(pm::_PM.AbstractPowerModel, x_list::Array{JuMP.VariableRef}, points; nw=0, id=1, var_name="x")

    gen_lines = []
    for j in 2:length(points)
        x1 = points[j-1].mw
        y1 = points[j-1].cost
        x2 = points[j].mw
        y2 = points[j].cost

        m = (y2 - y1)/(x2 - x1)

        if !isnan(m)
            b = y1 - m * x1
        else
            @assert isapprox(y1, y2)
            m = 0.0
            b = y1
        end
        push!(gen_lines, (slope=m, intercept=b))
    end

    pg_cost = JuMP.@variable(pm.model,
        base_name="$(nw)_pg_$(id)_cost",
        lower_bound = points[1].cost,
        upper_bound = points[end].cost
    )

    for line in gen_lines
        JuMP.@constraint(pm.model, pg_cost >= line.slope*var(pm, nw, :pg, id) + line.intercept)
    end

    return pg_cost

end

