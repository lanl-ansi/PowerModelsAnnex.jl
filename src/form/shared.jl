

NLForms = Union{NLSOCWROAForm, NLACRForm}

""
function PMs.objective_min_polynomial_fuel_cost(pm::GenericPowerModel{T}) where T <: NLForms
    PMs.check_polynomial_cost_models(pm)

    from_idx = Dict()
    for (n, nw_ref) in nws(pm)
        from_idx[n] = Dict(arc[1] => arc for arc in nw_ref[:arcs_from_dc])
    end

    pg = Dict()
    for n in nw_ids(pm), c in conductor_ids(pm, n), i in ids(pm, n, :gen)
        pg[(n,c,i)] = var(pm, n, c, :pg, i)
    end

    p_dc = Dict()
    for n in nw_ids(pm), c in conductor_ids(pm, n), i in ids(pm, n, :dcline)
        p_dc[(n,c,i)] = var(pm, n, c, :p_dc, from_idx[n][i])
    end

    return @NLobjective(pm.model, Min,
        sum(
            sum(   gen["cost"][1]*sum( pg[(n,c,i)] for c in conductor_ids(pm, n))^2 +
                   gen["cost"][2]*sum( pg[(n,c,i)] for c in conductor_ids(pm, n))+
                   gen["cost"][3] for (i,gen) in nw_ref[:gen]) +
            sum(   dcline["cost"][1]*sum( p_dc[(n,c,i)] for c in conductor_ids(pm, n))^2 +
                   dcline["cost"][2]*sum( p_dc[(n,c,i)] for c in conductor_ids(pm, n)) +
                   dcline["cost"][3] for (i,dcline) in nw_ref[:dcline])
        for (n, nw_ref) in nws(pm))
    )
end