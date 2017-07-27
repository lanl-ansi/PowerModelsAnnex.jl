

# Defines a variant of the SOCWRForm which is appropriate for outer approximation

export SOCWROAPowerModel, SOCWROAForm

@compat abstract type SOCWROAForm <: PMs.SOCWRForm end

const SOCWROAPowerModel = PMs.GenericPowerModel{SOCWROAForm}

"default SOCWROA constructor"
SOCWROAPowerModel(data::Dict{String,Any}; kwargs...) = PMs.GenericPowerModel(data, SOCWROAForm; kwargs...)

function PMs.constraint_voltage{T <: SOCWROAForm}(pm::GenericPowerModel{T})
    w = getindex(pm.model, :w)
    wr = getindex(pm.model, :wr)
    wi = getindex(pm.model, :wi)

    for (i,j) in keys(pm.ref[:buspairs])
        @NLconstraint(pm.model, (wr[(i,j)]^2 + wi[(i,j)]^2)/w[j] <= w[i])
    end
end

function PMs.constraint_thermal_limit_from{T <: SOCWROAForm}(pm::GenericPowerModel{T}, f_idx, rate_a)
    p_fr = getindex(pm.model, :p)[f_idx]
    q_fr = getindex(pm.model, :q)[f_idx]
    c = @NLconstraint(pm.model, sqrt(p_fr^2 + q_fr^2) <= rate_a)
    return Set([c])
end

function PMs.constraint_thermal_limit_to{T <: SOCWROAForm}(pm::GenericPowerModel{T}, t_idx, rate_a)
    p_to = getindex(pm.model, :p)[t_idx]
    q_to = getindex(pm.model, :q)[t_idx]
    c = @NLconstraint(pm.model, sqrt(p_to^2 + q_to^2) <= rate_a)
    return Set([c])
end

