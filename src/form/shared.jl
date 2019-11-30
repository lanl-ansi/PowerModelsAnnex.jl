NLPowerModels = Union{NLSOCWRPowerModel, NLACRPowerModel}

""
function PMs._objective_min_fuel_and_flow_cost_polynomial_linquad(pm::NLPowerModels)
    gen_cost = Dict()
    dcline_cost = Dict()

    for (n, nw_ref) in PMs.nws(pm)
        for (i,gen) in nw_ref[:gen]
            pg = sum( var(pm, n, c, :pg, i) for c in PMs.conductor_ids(pm, n) )

            if length(gen["cost"]) == 1
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1])
            elseif length(gen["cost"]) == 2
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*pg + gen["cost"][2])
            elseif length(gen["cost"]) == 3
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*pg^2 + gen["cost"][2]*pg + gen["cost"][3])
            else
                gen_cost[(n,i)] = 0.0
            end
        end

        from_idx = Dict(arc[1] => arc for arc in nw_ref[:arcs_from_dc])
        for (i,dcline) in nw_ref[:dcline]
            p_dc = sum( var(pm, n, c, :p_dc, from_idx[i]) for c in PMs.conductor_ids(pm, n) )

            if length(dcline["cost"]) == 1
                dcline_cost[(n,i)] = JuMP.@NLexpression(pm.model, dcline["cost"][1])
            elseif length(dcline["cost"]) == 2
                dcline_cost[(n,i)] = JuMP.@NLexpression(pm.model, dcline["cost"][1]*p_dc + dcline["cost"][2])
            elseif length(dcline["cost"]) == 3
                dcline_cost[(n,i)] = JuMP.@NLexpression(pm.model, dcline["cost"][1]*p_dc^2 + dcline["cost"][2]*p_dc + dcline["cost"][3])
            else
                dcline_cost[(n,i)] = 0.0
            end
        end
    end

    return JuMP.@NLobjective(pm.model, Min,
        sum(
            sum(    gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] ) +
            sum( dcline_cost[(n,i)] for (i,dcline) in nw_ref[:dcline] )
        for (n, nw_ref) in PMs.nws(pm))
    )
end
