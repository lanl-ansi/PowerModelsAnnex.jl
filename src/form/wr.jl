

# Defines a variant of the SOCWRForm which is appropriate for outer approximation

export SOCWROAPowerModel, SOCWROAForm

@compat abstract type SOCWROAForm <: PMs.SOCWRForm end

const SOCWROAPowerModel = PMs.GenericPowerModel{SOCWROAForm}

"default SOCWROA constructor"
SOCWROAPowerModel(data::Dict{String,Any}; kwargs...) = PMs.GenericPowerModel(data, SOCWROAForm; kwargs...)

""
function PMs.objective_min_fuel_cost{T <: SOCWROAForm}(pm::GenericPowerModel{T})
    PMs.check_cost_models(pm)

    pg = pm.var[:pg]
    dc_p = pm.var[:p_dc]

    from_idx = Dict(arc[1] => arc for arc in pm.ref[:arcs_from_dc])

    pg_sqr = pm.var[:pg_sqr] = @variable(pm.model, 
        [i in keys(pm.ref[:gen])], basename="pg_sqr",
        lowerbound = pm.ref[:gen][i]["pmin"]^2,
        upperbound = pm.ref[:gen][i]["pmax"]^2
    )
    for (i, gen) in pm.ref[:gen]
        @NLconstraint(pm.model, sqrt((2*pg[i])^2 + (pg_sqr[i]-1)^2) <= pg_sqr[i]+1)
    end

    dc_p_sqr = pm.var[:dc_p_sqr] = @variable(pm.model, 
        dc_p_sqr[i in keys(pm.ref[:dcline])], basename="dc_p_sqr",
        lowerbound = pm.ref[:dcline][i]["pminf"]^2,
        upperbound = pm.ref[:dcline][i]["pmaxf"]^2
    )
    for (i, dcline) in pm.ref[:dcline]
        @NLconstraint(pm.model, sqrt((2*dc_p[from_idx[i]])^2 + (dc_p_sqr[i]-1)^2) <= dc_p_sqr[i]+1)
    end

    return @objective(pm.model, Min,
        sum( gen["cost"][1]*pg_sqr[i] + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in pm.ref[:gen]) +
        sum(dcline["cost"][1]*dc_p_sqr[i]^2 + dcline["cost"][2]*dc_p[from_idx[i]] + dcline["cost"][3] for (i,dcline) in pm.ref[:dcline])
    )
end

function PMs.constraint_voltage{T <: SOCWROAForm}(pm::GenericPowerModel{T})
    w = pm.var[:w]
    wr = pm.var[:wr]
    wi = pm.var[:wi]

    for (i,j) in keys(pm.ref[:buspairs])
        @NLconstraint(pm.model, (wr[(i,j)]^2 + wi[(i,j)]^2)/w[j] <= w[i])
    end
end

function PMs.constraint_thermal_limit_from{T <: SOCWROAForm}(pm::GenericPowerModel{T}, f_idx, rate_a)
    p_fr = pm.var[:p][f_idx]
    q_fr = pm.var[:q][f_idx]
    @NLconstraint(pm.model, sqrt(p_fr^2 + q_fr^2) <= rate_a)
end

function PMs.constraint_thermal_limit_to{T <: SOCWROAForm}(pm::GenericPowerModel{T}, t_idx, rate_a)
    p_to = pm.var[:p][t_idx]
    q_to = pm.var[:q][t_idx]
    @NLconstraint(pm.model, sqrt(p_to^2 + q_to^2) <= rate_a)
end

