export run_opf_pwl_delta, build_opf_pwl_delta

"computes the cost values that are relivent, based on min/max values"
function get_active_cost_points(component)
    pmin = component["pmin"]
    pmax = component["pmax"]
    ncost = component["ncost"]
    points = component["cost"]

    first_active = 1
    for i in 1:(ncost-1)
        x0 = points[2*i-1]
        x1 = points[2*(i+1)-1]
        if pmin >= x0
            first_active = i
        end
        if pmin <= x1
            break
        end
    end

    last_active = ncost
    for i in 1:(ncost-1)
        x0 = points[end - (2*(i+1)-1)]
        x1 = points[end - (2*i-1)]
        if pmax <= x1
            last_active = ncost - i + 1
        end
        if pmax >= x0
            break
        end
    end

    points = points[2*first_active - 1 : 2*last_active]
    ncost = div(length(points), 2)

    @assert points[1] <= pmin && points[1+2] >= pmin
    @assert points[end-3] <= pmax && points[end-1] >= pmax

    return ncost, points
end


function calc_comp_lines_from_points(points)
    line_data = []
    for i in 3:2:length(points)
        x1 = points[i-2]
        y1 = points[i-1]
        x2 = points[i-0]
        y2 = points[i+1]

        m = (y2 - y1)/(x2 - x1)
        b = y1 - m * x1

        push!(line_data, (slope=m, intercept=b))
    end

    for i in 2:length(line_data)
        if line_data[i-1].slope > line_data[i].slope
            Memento.error(_LOGGER, "non-convex pwl function found in points $(points)\nlines: $(line_data)")
        end
    end

    return line_data
end


""
function run_opf_pwl_delta(file, model_type::Type, optimizer; kwargs...)
    return _PM.run_model(file, model_type, optimizer, build_opf_pwl_delta; kwargs...)
end

"a variant of the OPF problem specification that uses a max variant of the pwl cost function implementation"
function build_opf_pwl_delta(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_delta(pm)

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
function objective_min_fuel_and_flow_cost_delta(pm::_PM.AbstractPowerModel; kwargs...)
    model = _PM.check_cost_models(pm)

    if model == 1
        return objective_min_fuel_and_flow_cost_pwl_delta(pm; kwargs...)
    elseif model == 2
        return objective_min_fuel_and_flow_cost_polynomial(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end


""
function objective_min_fuel_and_flow_cost_pwl_delta(pm::_PM.AbstractPowerModel; kwargs...)
    objective_variable_pg_cost_delta(pm; kwargs...)
    _PM.objective_variable_dc_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_cost, i) for (i,gen) in nw_ref[:gen]) +
            sum( var(pm, n, :p_dc_cost, i) for (i,dcline) in nw_ref[:dcline])
        for (n, nw_ref) in _PM.nws(pm))
    )
end


"adds pg_cost variables and constraints"
function objective_variable_pg_cost_delta(pm::_PM.AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in _PM.nws(pm)
        pg_cost = var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in ref(pm, n, :gen)
            pmin = gen["pmin"]
            mws = Float64[]
            costs = Float64[]

            ncost, points = get_active_cost_points(gen)
            for i in 1:ncost
                push!(mws, points[2*i-1])
                push!(costs, points[2*i])
            end

            cost_per_mw = Float64[0.0]
            for i in 2:ncost
                x0 = mws[i-1]
                y0 = costs[i-1]
                x1 = mws[i]
                y1 = costs[i]

                m = (y1 - y0)/(x1 - x0)
                if !isnan(m)
                    push!(cost_per_mw, m)
                else
                    @assert isapprox(y0, y1)
                    push!(cost_per_mw, 0.0)
                end
            end

            pg_cost_mw = JuMP.@variable(pm.model,
                [i in 2:ncost], base_name="$(n)_pg_cost_mw",
                lower_bound = 0.0,
                upper_bound = mws[i] - mws[i-1]
            )

            JuMP.@constraint(pm.model, mws[1] + sum(pg_cost_mw[i] for i in 2:ncost) == sum(var(pm, n, :pg, i)[c] for c in _PM.conductor_ids(pm, n)))
            pg_cost[i] = costs[1] + sum(cost_per_mw[i]*pg_cost_mw[i] for i in 2:ncost)
        end

        report && _IM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)
    end
end
