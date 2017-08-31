

# Defines a variant of the SOCWRForm which is appropriate for outer approximation

export SOCWROAPowerModel, SOCWROAForm

@compat abstract type SOCWROAForm <: PMs.SOCWRForm end

const SOCWROAPowerModel = PMs.GenericPowerModel{SOCWROAForm}

"default SOCWROA constructor"
SOCWROAPowerModel(data::Dict{String,Any}; kwargs...) = PMs.GenericPowerModel(data, SOCWROAForm; kwargs...)

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

