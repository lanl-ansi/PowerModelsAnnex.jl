

# Defines a variant of the SOCWRForm which is appropriate for outer approximation

export SOCWROAPowerModel, SOCWROAForm

abstract type SOCWROAForm <: PMs.SOCWRForm end

const SOCWROAPowerModel = PMs.GenericPowerModel{SOCWROAForm}

"default SOCWROA constructor"
SOCWROAPowerModel(data::Dict{String,Any}; kwargs...) = PMs.GenericPowerModel(data, SOCWROAForm; kwargs...)

""
function PMs._objective_min_fuel_and_flow_cost_polynomial_linquad(pm::PMs.GenericPowerModel{T}) where T <: SOCWROAForm
    gen_cost = Dict()
    dcline_cost = Dict()

    for (n, nw_ref) in PMs.nws(pm)

        var(pm, n)[:pg_sqr] = Dict()
        for (i,gen) in nw_ref[:gen]
            pg = sum( var(pm, n, c, :pg, i) for c in PMs.conductor_ids(pm, n) )

            if length(gen["cost"]) == 1
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1])
            elseif length(gen["cost"]) == 2
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*pg + gen["cost"][2])
            elseif length(gen["cost"]) == 3
                pmin = sum(gen["pmin"][c] for c in PMs.conductor_ids(pm, n))
                pmax = sum(gen["pmax"][c] for c in PMs.conductor_ids(pm, n))

                pg_sqr_ub = max(pmin^2, pmax^2)
                pg_sqr_lb = 0.0
                if pmin > 0.0
                    pg_sqr_lb = pmin^2
                end
                if pmax < 0.0
                    pg_sqr_lb = pmax^2
                end

                pg_sqr = var(pm, n, :pg_sqr)[i] = JuMP.@variable(pm.model,
                    base_name="$(n)_pg_sqr_$(i)",
                    lower_bound = pg_sqr_lb,
                    upper_bound = pg_sqr_ub,
                    start = 0.0
                )
                JuMP.@NLconstraint(pm.model, pg^2 <= pg_sqr)

                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*pg_sqr + gen["cost"][2]*pg + gen["cost"][3])
            else
                gen_cost[(n,i)] = 0.0
            end
        end

        from_idx = Dict(arc[1] => arc for arc in nw_ref[:arcs_from_dc])

        var(pm, n)[:p_dc_sqr] = Dict()
        for (i,dcline) in nw_ref[:dcline]
            p_dc = sum( var(pm, n, c, :p_dc, from_idx[i]) for c in PMs.conductor_ids(pm, n) )

            if length(dcline["cost"]) == 1
                dcline_cost[(n,i)] = JuMP.@NLexpression(pm.model, dcline["cost"][1])
            elseif length(dcline["cost"]) == 2
                dcline_cost[(n,i)] = JuMP.@NLexpression(pm.model, dcline["cost"][1]*p_dc + dcline["cost"][2])
            elseif length(dcline["cost"]) == 3
                pmin = sum(dcline["pminf"][c] for c in PMs.conductor_ids(pm, n))
                pmax = sum(dcline["pmaxf"][c] for c in PMs.conductor_ids(pm, n))

                p_dc_sqr_ub = max(pmin^2, pmax^2)
                p_dc_sqr_lb = 0.0
                if pmin > 0.0
                    p_dc_sqr_lb = pmin^2
                end
                if pmax < 0.0
                    p_dc_sqr_lb = pmax^2
                end

                p_dc_sqr = var(pm, n, :p_dc_sqr)[i] = JuMP.@variable(pm.model,
                    base_name="$(n)_p_dc_sqr_$(i)",
                    lower_bound = p_dc_sqr_lb,
                    upper_bound = p_dc_sqr_ub,
                    start = 0.0
                )
                JuMP.@NLconstraint(pm.model, p_dc^2 <= p_dc_sqr)

                dcline_cost[(n,i)] = JuMP.@NLexpression(pm.model, dcline["cost"][1]*p_dc_sqr + dcline["cost"][2]*p_dc + dcline["cost"][3])
            else
                dcline_cost[(n,i)] = 0.0
            end
        end
    end

    return JuMP.@NLobjective(pm.model, Min,
        sum(
            sum( gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] ) +
            sum( dcline_cost[(n,i)] for (i,dcline) in nw_ref[:dcline] )
        for (n, nw_ref) in PMs.nws(pm))
    )
end


function PMs.constraint_model_voltage(pm::PMs.GenericPowerModel{T}, n::Int, h::Int) where T <: SOCWROAForm
    w  = var(pm, n, h,  :w)
    wr = var(pm, n, h, :wr)
    wi = var(pm, n, h, :wi)


    for (i,j) in PMs.ids(pm, n, :buspairs)
        @NLconstraint(pm.model, (wr[(i,j)]^2 + wi[(i,j)]^2)/w[j] <= w[i])
    end
end

function PMs.constraint_thermal_limit_from(pm::PMs.GenericPowerModel{T}, n::Int, h::Int, f_idx, rate_a) where T <: SOCWROAForm
    p_fr = var(pm, n, h, :p, f_idx)
    q_fr = var(pm, n, h, :q, f_idx)
    @NLconstraint(pm.model, sqrt(p_fr^2 + q_fr^2) <= rate_a)
end

function PMs.constraint_thermal_limit_to(pm::PMs.GenericPowerModel{T}, n::Int, h::Int, t_idx, rate_a) where T <: SOCWROAForm
    p_to = var(pm, n, h, :p, t_idx)
    q_to = var(pm, n, h, :q, t_idx)
    @NLconstraint(pm.model, sqrt(p_to^2 + q_to^2) <= rate_a)
end



# Defines a variant of the SOCWRForm which forces the NL solver path

export NLSOCWRPowerModel, NLSOCWROAForm

abstract type NLSOCWROAForm <: PMs.SOCWRForm end

const NLSOCWRPowerModel = PMs.GenericPowerModel{NLSOCWROAForm}

"default NLSOCWROAForm constructor"
NLSOCWRPowerModel(data::Dict{String,Any}; kwargs...) = PMs.GenericPowerModel(data, NLSOCWROAForm; kwargs...)



# Defines a variant of the QCWRTriForm without the linking constraints

export QCWRTriNoLinkPowerModel, QCWRTriNoLinkForm

abstract type QCWRTriNoLinkForm <: PMs.QCWRTriForm end

const QCWRTriNoLinkPowerModel = PMs.GenericPowerModel{QCWRTriNoLinkForm}

"default QC trilinear without linking constraint model constructor"
QCWRTriNoLinkPowerModel(data::Dict{String,Any}; kwargs...) = PMs.GenericPowerModel(data, QCWRTriNoLinkForm; kwargs...)

function PMs.constraint_model_voltage(pm::PMs.GenericPowerModel{T}, n::Int, c::Int) where T <: QCWRTriNoLinkForm
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

    for bp in PMs.ids(pm, n, :buspairs)
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
