NLPowerModels = Union{NLSOCWRPowerModel, NLACRPowerModel}

"force JuMP to define an NLobjective"
function _PM.objective_min_fuel_and_flow_cost(pm::NLPowerModels; kwargs...)
    _PM.expression_pg_cost(pm; kwargs...)
    _PM.expression_p_dc_cost(pm; kwargs...)

    return JuMP.@objective(pm.model, Min,
        sum(
            sum( var(pm, n,   :pg_cost, i) for (i,gen) in nw_ref[:gen]) +
            sum( var(pm, n, :p_dc_cost, i) for (i,dcline) in nw_ref[:dcline])
        for (n, nw_ref) in _PM.nws(pm))
    )
end
