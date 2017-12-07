

# Defines a variant of the SOCWRForm which is appropriate for outer approximation

export SOCWROAPowerModel, SOCWROAForm

@compat abstract type SOCWROAForm <: PMs.SOCWRForm end

const SOCWROAPowerModel = PMs.GenericPowerModel{SOCWROAForm}

"default SOCWROA constructor"
SOCWROAPowerModel(data::Dict{String,Any}; kwargs...) = PMs.GenericPowerModel(data, SOCWROAForm; kwargs...)

""
function PMs.objective_min_fuel_cost{T <: SOCWROAForm}(pm::GenericPowerModel{T}, nws=[pm.cnw])
    @assert !(pm.data["multinetwork"])

    PMs.check_cost_models(pm, nws)

    pg = var(pm, :pg)
    dc_p = var(pm, :p_dc)

    from_idx = Dict(arc[1] => arc for arc in ref(pm, :arcs_from_dc))

    pm.var[:pg_sqr] = Dict{Int, Any}()
    @expression(pm.model, gen_cost, 0)
    for (i, gen) in ref(pm, :gen)
        if gen["cost"][1] != 0.0
            pg_sqr = pm.var[:pg_sqr][i] = @variable(pm.model, 
                basename="pg_sqr",
                lowerbound = ref(pm, :gen)[i]["pmin"]^2,
                upperbound = ref(pm, :gen)[i]["pmax"]^2
            )
            @NLconstraint(pm.model, sqrt((2*pg[i])^2 + (pg_sqr-1)^2) <= pg_sqr+1)
            gen_cost = gen_cost + gen["cost"][1]*pg_sqr + gen["cost"][2]*pg[i] + gen["cost"][3]
        else
            gen_cost = gen_cost + gen["cost"][2]*pg[i] + gen["cost"][3]
        end
    end

    pm.var[:dc_p_sqr] = Dict{Int, Any}()
    @expression(pm.model, dcline_cost, 0)
    for (i, dcline) in ref(pm, :dcline)
        if dcline["cost"][1] != 0.0
            dc_p_sqr = pm.var[:dc_p_sqr][i] = @variable(pm.model, 
                basename="dc_p_sqr",
                lowerbound = ref(pm, :dcline)[i]["pminf"]^2,
                upperbound = ref(pm, :dcline)[i]["pmaxf"]^2
            )
            @NLconstraint(pm.model, sqrt((2*dc_p[from_idx[i]])^2 + (dc_p_sqr-1)^2) <= dc_p_sqr+1)
            dcline_cost = dcline_cost + dcline["cost"][1]*dc_p_sqr^2 + dcline["cost"][2]*dc_p[from_idx[i]] + dcline["cost"][3]
        else
            dcline_cost = dcline_cost + dcline["cost"][2]*dc_p[from_idx[i]] + dcline["cost"][3]
        end
    end

    return @objective(pm.model, Min, gen_cost + dcline_cost)
end


function PMs.constraint_voltage{T <: SOCWROAForm}(pm::GenericPowerModel{T}, n::Int)
    w = pm.var[:nw][n][:w]
    wr = pm.var[:nw][n][:wr]
    wi = pm.var[:nw][n][:wi]

    for (i,j) in keys(pm.ref[:nw][n][:buspairs])
        @NLconstraint(pm.model, (wr[(i,j)]^2 + wi[(i,j)]^2)/w[j] <= w[i])
    end
end

function PMs.constraint_thermal_limit_from{T <: SOCWROAForm}(pm::GenericPowerModel{T}, n::Int, f_idx, rate_a)
    p_fr = pm.var[:nw][n][:p][f_idx]
    q_fr = pm.var[:nw][n][:q][f_idx]
    @NLconstraint(pm.model, sqrt(p_fr^2 + q_fr^2) <= rate_a)
end

function PMs.constraint_thermal_limit_to{T <: SOCWROAForm}(pm::GenericPowerModel{T}, n::Int, t_idx, rate_a)
    p_to = pm.var[:nw][n][:p][t_idx]
    q_to = pm.var[:nw][n][:q][t_idx]
    @NLconstraint(pm.model, sqrt(p_to^2 + q_to^2) <= rate_a)
end



# Defines a variant of the SOCWRForm which forces the NL solver path

export NLSOCWRPowerModel, NLSOCWROAForm

@compat abstract type NLSOCWROAForm <: PMs.SOCWRForm end

const NLSOCWRPowerModel = PMs.GenericPowerModel{NLSOCWROAForm}

"default NLSOCWROAForm constructor"
NLSOCWRPowerModel(data::Dict{String,Any}; kwargs...) = PMs.GenericPowerModel(data, NLSOCWROAForm; kwargs...)





