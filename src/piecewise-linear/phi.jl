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

    objective_min_fuel_and_flow_cost_phi(pm)

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
function objective_min_fuel_and_flow_cost_phi(pm::_PM.AbstractPowerModel; kwargs...)
    model = _PM.check_cost_models(pm)

    if model == 1
        return objective_min_fuel_and_flow_cost_pwl_phi(pm; kwargs...)
    elseif model == 2
        return objective_min_fuel_and_flow_cost_polynomial(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end
end

""
function objective_min_fuel_and_flow_cost_pwl_phi(pm::_PM.AbstractPowerModel; kwargs...)
    objective_variable_pg_cost_phi(pm; kwargs...)
    _PM.objective_variable_dc_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_cost, i) for (i,gen) in nw_ref[:gen]) +
            sum( var(pm, n, :p_dc_cost, i) for (i,dcline) in nw_ref[:dcline])
        for (n, nw_ref) in _PM.nws(pm))
    )
end

"adds pg_cost variables and constraints"
function objective_variable_pg_cost_phi(pm::_PM.AbstractPowerModel, report::Bool=true)
    for (n, nw_ref) in _PM.nws(pm)
        pg_cost = var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in ref(pm, n, :gen)
            ncost, points = get_active_cost_points(gen)

            mws = Float64[]
            costs = Float64[]

            for j in 1:ncost
                push!(mws, points[2*j-1])
                push!(costs, points[2*j])
            end


            cost_per_mw = Float64[0.0]
            cost_per_mw_b = Float64[0.0]
            for j in 2:ncost
                x0 = mws[j-1]
                y0 = costs[j-1]
                x1 = mws[j]
                y1 = costs[j]

                m = (y1 - y0)/(x1 - x0)
                if !isnan(m)
                    b = y1 - m * x1

                    push!(cost_per_mw, m)
                    push!(cost_per_mw_b, b)
                else
                    #println(y0, " ", y1)
                    @assert isapprox(y0, y1)
                    push!(cost_per_mw, 0.0)
                    push!(cost_per_mw_b, y0)
                end
            end

            # println(cost_per_mw)
            # println(mws)
            # println()

            pmax = gen["pmax"]
            pg_phi = JuMP.@variable(pm.model,
                [j in 3:ncost], base_name="$(n)_pg_phi",
                lower_bound = 0.0,
                upper_bound = pmax - mws[j-1]
            )

            for j in 3:ncost
                JuMP.@constraint(pm.model, pg_phi[j] >= sum(var(pm, n, :pg, i)[c] for c in _PM.conductor_ids(pm, n)) - mws[j-1])
            end

            #pg_cost[i] = cost_per_mw[2] * (sum(var(pm, n, :pg, i)[c] for c in _PM.conductor_ids(pm, n)) - mws[1]) + costs[1]
            pg_cost[i] = cost_per_mw[2] * (sum(var(pm, n, :pg, i)[c] for c in _PM.conductor_ids(pm, n))) + cost_per_mw_b[2]
            if ncost > 2
                pg_cost[i] += sum((cost_per_mw[j] - cost_per_mw[j-1])*pg_phi[j] for j in 3:ncost)
            end
        end

        report && _IM.sol_component_value(pm, n, :gen, :pg_cost, ids(pm, n, :gen), pg_cost)
    end
end
