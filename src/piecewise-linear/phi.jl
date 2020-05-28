export run_opf_pwl_phi, build_opf_pwl_phi

""
function run_opf_pwl_phi(file, model_type::Type, optimizer; kwargs...)
    return _PM.run_model(file, model_type, optimizer, build_opf_pwl_phi; kwargs...)
end

"a variant of the OPF problem specification that uses a max variant of the pwl cost function implementation"
function build_opf_pwl_phi(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_dcline_power(pm)

    objective_min_fuel_cost_phi(pm)

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
function objective_min_fuel_cost_phi(pm::_PM.AbstractPowerModel; kwargs...)
    model = _PM.check_cost_models(pm)

    if model == 1
        return objective_min_fuel_cost_pwl_phi(pm; kwargs...)
    elseif model == 2
        return _PM.objective_min_fuel_cost_polynomial(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end
end

""
function objective_min_fuel_cost_pwl_phi(pm::_PM.AbstractPowerModel; kwargs...)
    objective_variable_pg_cost_phi(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_cost, i) for (i,gen) in nw_ref[:gen])
        for (n, nw_ref) in _PM.nws(pm))
    )
end

"adds pg_cost variables and constraints"
function objective_variable_pg_cost_phi(pm::_PM.AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in _PM.nws(pm)
        pg_cost = var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in ref(pm, n, :gen)
            points = _PM.calc_pwl_points(gen["ncost"], gen["cost"], gen["pmin"], gen["pmax"])

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
                    @assert isapprox(y0, y1)
                    m = 0.0
                    b = y0
                end
                push!(gen_lines, (slope=m, intercept=b))
            end

            pmax = gen["pmax"]
            pg_phi = JuMP.@variable(pm.model,
                [j in 3:length(points)], base_name="$(n)_pg_phi",
                lower_bound = 0.0,
                upper_bound = pmax - points[j-1].mw
            )

            for j in 3:length(points)
                JuMP.@constraint(pm.model, pg_phi[j] >= sum(var(pm, n, :pg, i)[c] for c in _PM.conductor_ids(pm, n)) - points[j-1].mw)
            end

            pg_cost[i] = gen_lines[1].slope * (sum(var(pm, n, :pg, i)[c] for c in _PM.conductor_ids(pm, n))) + gen_lines[1].intercept
            if length(points) > 2
                pg_cost[i] += sum((gen_lines[j-1].slope - gen_lines[j-2].slope)*pg_phi[j] for j in 3:length(points))
            end
        end

        report && _IM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)
    end
end
