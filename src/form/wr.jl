

# Defines a variant of the SOCWRForm which is appropriate for outer approximation

export SOCWROAPowerModel, SOCWROAForm

abstract type SOCWROAForm <: PMs.SOCWRForm end

const SOCWROAPowerModel = PMs.GenericPowerModel{SOCWROAForm}

"default SOCWROA constructor"
SOCWROAPowerModel(data::Dict{String,Any}; kwargs...) = PMs.GenericPowerModel(data, SOCWROAForm; kwargs...)

""
function PMs.objective_min_fuel_cost{T <: SOCWROAForm}(pm::GenericPowerModel{T})
    @assert !InfrastructureModels.ismultinetwork(pm.data)
    @assert !haskey(pm.data, "conductors")

    PMs.check_cost_models(pm)

    pg = var(pm, :pg)
    dc_p = var(pm, :p_dc)

    from_idx = Dict(arc[1] => arc for arc in ref(pm, :arcs_from_dc))

    pm.var[:pg_sqr] = Dict{Int, Any}()
    @expression(pm.model, gen_cost, 0)
    for (i, gen) in ref(pm, :gen)
        if gen["cost"][1] != 0.0
            pg_sqr = pm.var[:pg_sqr][i] = @variable(pm.model, 
                basename="pg_sqr",
                lowerbound = ref(pm, :gen, i, "pmin")^2,
                upperbound = ref(pm, :gen, i, "pmax")^2
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
                lowerbound = ref(pm, :dcline, i, "pminf")^2,
                upperbound = ref(pm, :dcline, i, "pmaxf")^2
            )
            @NLconstraint(pm.model, sqrt((2*dc_p[from_idx[i]])^2 + (dc_p_sqr-1)^2) <= dc_p_sqr+1)
            dcline_cost = dcline_cost + dcline["cost"][1]*dc_p_sqr^2 + dcline["cost"][2]*dc_p[from_idx[i]] + dcline["cost"][3]
        else
            dcline_cost = dcline_cost + dcline["cost"][2]*dc_p[from_idx[i]] + dcline["cost"][3]
        end
    end

    return @objective(pm.model, Min, gen_cost + dcline_cost)
end


function PMs.constraint_voltage(pm::GenericPowerModel{T}, n::Int, h::Int) where T <: SOCWROAForm
    w  = var(pm, n, h,  :w)
    wr = var(pm, n, h, :wr)
    wi = var(pm, n, h, :wi)


    for (i,j) in ids(pm, n, :buspairs)
        @NLconstraint(pm.model, (wr[(i,j)]^2 + wi[(i,j)]^2)/w[j] <= w[i])
    end
end

function PMs.constraint_thermal_limit_from{T <: SOCWROAForm}(pm::GenericPowerModel{T}, n::Int, h::Int, f_idx, rate_a)
    p_fr = var(pm, n, h, :p, f_idx)
    q_fr = var(pm, n, h, :q, f_idx)
    @NLconstraint(pm.model, sqrt(p_fr^2 + q_fr^2) <= rate_a)
end

function PMs.constraint_thermal_limit_to{T <: SOCWROAForm}(pm::GenericPowerModel{T}, n::Int, h::Int, t_idx, rate_a)
    p_to = var(pm, n, h, :p, t_idx)
    q_to = var(pm, n, h, :q, t_idx)
    @NLconstraint(pm.model, sqrt(p_to^2 + q_to^2) <= rate_a)
end



# Defines a variant of the SOCWRForm which forces the NL solver path

export NLSOCWRPowerModel, NLSOCWROAForm

@compat abstract type NLSOCWROAForm <: PMs.SOCWRForm end

const NLSOCWRPowerModel = PMs.GenericPowerModel{NLSOCWROAForm}

"default NLSOCWROAForm constructor"
NLSOCWRPowerModel(data::Dict{String,Any}; kwargs...) = PMs.GenericPowerModel(data, NLSOCWROAForm; kwargs...)



# Defines a variant of the QCWRTriForm without the linking constraints

export QCWRTriNoLinkPowerModel, QCWRTriNoLinkForm

@compat abstract type QCWRTriNoLinkForm <: PMs.QCWRTriForm end

const QCWRTriNoLinkPowerModel = PMs.GenericPowerModel{QCWRTriNoLinkForm}

"default QC trilinear without linking constraint model constructor"
QCWRTriNoLinkPowerModel(data::Dict{String,Any}; kwargs...) = PMs.GenericPowerModel(data, QCWRTriNoLinkForm; kwargs...)

function PMs.constraint_voltage(pm::GenericPowerModel{T}, n::Int, c::Int) where T <: QCWRTriNoLinkForm
    v = var(pm, n, c, :vm)
    t = var(pm, n, c, :va)

    td = var(pm, n, c, :td)
    si = var(pm, n, c, :si)
    cs = var(pm, n, c, :cs)

    w = var(pm, n, c, :w)
    wr = var(pm, n, c, :wr)
    lambda_wr = var(pm, n, c, :lambda_wr)
    wi = var(pm, n, c, :wi)
    lambda_wi = var(pm, n, c, :lambda_wi)

    for (i,b) in ref(pm, n, :bus)
        InfrastructureModels.relaxation_sqr(pm.model, v[i], w[i])
    end

    for bp in ids(pm, n, :buspairs)
        i,j = bp
        @constraint(pm.model, t[i] - t[j] == td[bp])

        PMs.relaxation_sin(pm.model, td[bp], si[bp])
        PMs.relaxation_cos(pm.model, td[bp], cs[bp])
        InfrastructureModels.relaxation_trilinear(pm.model, v[i], v[j], cs[bp], wr[bp], lambda_wr[bp,:])
        InfrastructureModels.relaxation_trilinear(pm.model, v[i], v[j], si[bp], wi[bp], lambda_wi[bp,:])

        # this constraint is redudant and useful for debugging
        #InfrastructureModels.relaxation_complex_product(pm.model, w[i], w[j], wr[bp], wi[bp])
   end

   for (i,branch) in ref(pm, n, :branch)
        pair = (branch["f_bus"], branch["t_bus"])
        buspair = ref(pm, n, :buspairs, pair)

        # to prevent this constraint from being posted on multiple parallel branchs
        if buspair["branch"] == i
            PMs.constraint_power_magnitude_sqr(pm, i, nw=n, cnd=c)
            PMs.constraint_power_magnitude_link(pm, i, nw=n, cnd=c)
        end
    end

end

