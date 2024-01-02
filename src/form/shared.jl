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


# ""
# function _PM._objective_min_fuel_and_flow_cost_polynomial_linquad(pm::NLPowerModels, report::Bool=true)
#     gen_cost = Dict()
#     dcline_cost = Dict()

#     for (n, nw_ref) in _PM.nws(pm)
#         for (i,gen) in nw_ref[:gen]
#             pg = sum( var(pm, n, :pg, i)[c] for c in _PM.conductor_ids(pm, n) )

#             if length(gen["cost"]) == 1
#                 gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1])
#             elseif length(gen["cost"]) == 2
#                 gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*pg + gen["cost"][2])
#             elseif length(gen["cost"]) == 3
#                 gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*pg^2 + gen["cost"][2]*pg + gen["cost"][3])
#             else
#                 gen_cost[(n,i)] = 0.0
#             end
#         end

#         from_idx = Dict(arc[1] => arc for arc in nw_ref[:arcs_from_dc])
#         for (i,dcline) in nw_ref[:dcline]
#             p_dc = sum( var(pm, n, :p_dc, from_idx[i])[c] for c in _PM.conductor_ids(pm, n) )

#             if length(dcline["cost"]) == 1
#                 dcline_cost[(n,i)] = JuMP.@NLexpression(pm.model, dcline["cost"][1])
#             elseif length(dcline["cost"]) == 2
#                 dcline_cost[(n,i)] = JuMP.@NLexpression(pm.model, dcline["cost"][1]*p_dc + dcline["cost"][2])
#             elseif length(dcline["cost"]) == 3
#                 dcline_cost[(n,i)] = JuMP.@NLexpression(pm.model, dcline["cost"][1]*p_dc^2 + dcline["cost"][2]*p_dc + dcline["cost"][3])
#             else
#                 dcline_cost[(n,i)] = 0.0
#             end
#         end
#     end

#     return JuMP.@NLobjective(pm.model, Min,
#         sum(
#             sum(    gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] ) +
#             sum( dcline_cost[(n,i)] for (i,dcline) in nw_ref[:dcline] )
#         for (n, nw_ref) in _PM.nws(pm))
#     )
# end
