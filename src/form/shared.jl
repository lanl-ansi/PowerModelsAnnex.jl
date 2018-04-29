

NLForms = Union{NLSOCWROAForm, NLACRForm}

""
function objective_min_polynomial_fuel_cost(pm::GenericPowerModel{T}) where T <: NLForms
    PMs.check_polynomial_cost_models(pm)

    from_idx = Dict()
    for (n, nw_ref) in nws(pm)
        from_idx[n] = Dict(arc[1] => arc for arc in nw_ref[:arcs_from_dc])
    end

    return @NLobjective(pm.model, Min,
        sum(
            sum(
                sum(   getmpv(gen["cost"],h)[1]*var(pm, n, h, :pg, i)^2 + getmpv(gen["cost"],h)[2]*var(pm, n, h, :pg, i) + getmpv(gen["cost"],h)[3] for (i,gen) in nw_ref[:gen]) +
                sum(getmpv(dcline["cost"],h)[1]*var(pm, n, h, :p_dc, from_idx[n][i])^2 + getmpv(dcline["cost"],h)[2]*var(pm, n, h, :p_dc, from_idx[n][i]) + getmpv(dcline["cost"],h)[3] for (i,dcline) in nw_ref[:dcline])
            for h in phase_ids(pm, n))
        for (n, nw_ref) in nws(pm))
    )
end