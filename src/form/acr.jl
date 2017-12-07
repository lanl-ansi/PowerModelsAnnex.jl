export NLACRPowerModel, NLACRForm

@compat abstract type NLACRForm <: PMs.AbstractACRForm end

const NLACRPowerModel = PMs.GenericPowerModel{NLACRForm}

"default NLACRForm constructor"
NLACRPowerModel(data::Dict{String,Any}; kwargs...) = PMs.GenericPowerModel(data, NLACRForm; kwargs...)

"a copy of the standard objective_min_fuel_cost, but with an NL objective to force NL path"
function PMs.objective_min_fuel_cost{T <: NLACRForm}(pm::GenericPowerModel{T}, nws=[pm.cnw])
    PMs.check_cost_models(pm, nws)

    pg = Dict(n => pm.var[:nw][n][:pg] for n in nws)
    dc_p = Dict(n => pm.var[:nw][n][:p_dc] for n in nws)

    from_idx = Dict()
    for n in nws
        ref = pm.ref[:nw][n]
        from_idx[n] = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])
    end

    return @NLobjective(pm.model, Min, 
        sum(
            sum(gen["cost"][1]*pg[n][i]^2 + gen["cost"][2]*pg[n][i] + gen["cost"][3] for (i,gen) in pm.ref[:nw][n][:gen]) +
            sum(dcline["cost"][1]*dc_p[n][from_idx[n][i]]^2 + dcline["cost"][2]*dc_p[n][from_idx[n][i]] + dcline["cost"][3] for (i,dcline) in pm.ref[:nw][n][:dcline])
        for n in nws)
    )
end
