NLPowerModels = Union{NLSOCWRPowerModel, NLACRPowerModel}

"force JuMP to define an NLobjective"
function _PM.objective_min_fuel_and_flow_cost(pm::NLPowerModels; kwargs...)
    nl_gen = _PM.check_nl_gen_cost_models(pm)
    nl_dc = _PM.check_nl_dcline_cost_models(pm)

    _PM.expression_pg_cost(pm; kwargs...)
    _PM.expression_p_dc_cost(pm; kwargs...)

    pg_cost = Dict()
    p_dc_cost = Dict()
    for (n, nw_ref) in _PM.nws(pm)
        for (i,gen) in nw_ref[:gen]
            pg_cost[(n,i)] = var(pm, n,   :pg_cost, i)
        end
        for (i,dcline) in nw_ref[:dcline]
            p_dc_cost[(n,i)] = var(pm, n, :p_dc_cost, i)
        end
    end

    return JuMP.@NLobjective(pm.model, Min,
        sum(
            sum( pg_cost[n,i] for (i,gen) in nw_ref[:gen]) +
            sum( p_dc_cost[n,i] for (i,dcline) in nw_ref[:dcline])
        for (n, nw_ref) in _PM.nws(pm))
    )
end
