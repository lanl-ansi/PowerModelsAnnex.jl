

# Defines a variant of the SOCWRForm which is appropriate for outer approximation

export SOCWROAPowerModel

mutable struct SOCWROAPowerModel <: _PM.AbstractSOCWRModel _PM.@pm_fields end


""
function _PM._objective_min_fuel_and_flow_cost_polynomial_linquad(pm::SOCWROAPowerModel, report::Bool=true)
    gen_cost = Dict()
    dcline_cost = Dict()

    for (n, nw_ref) in _PM.nws(pm)

        var(pm, n)[:pg_sqr] = Dict()
        for (i,gen) in nw_ref[:gen]
            pg = sum( var(pm, n, :pg, i)[c] for c in _PM.conductor_ids(pm, n) )

            if length(gen["cost"]) == 1
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1])
            elseif length(gen["cost"]) == 2
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*pg + gen["cost"][2])
            elseif length(gen["cost"]) == 3
                pmin = sum(gen["pmin"][c] for c in _PM.conductor_ids(pm, n))
                pmax = sum(gen["pmax"][c] for c in _PM.conductor_ids(pm, n))

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

                if report
                    sol(pm, n, :gen, i)[:pg_sqr] = pg_sqr
                end

                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*pg_sqr + gen["cost"][2]*pg + gen["cost"][3])
            else
                gen_cost[(n,i)] = 0.0
            end
        end

        from_idx = Dict(arc[1] => arc for arc in nw_ref[:arcs_from_dc])

        var(pm, n)[:p_dc_sqr] = Dict()
        for (i,dcline) in nw_ref[:dcline]
            p_dc = sum( var(pm, n, :p_dc, from_idx[i])[c] for c in _PM.conductor_ids(pm, n) )

            if length(dcline["cost"]) == 1
                dcline_cost[(n,i)] = JuMP.@NLexpression(pm.model, dcline["cost"][1])
            elseif length(dcline["cost"]) == 2
                dcline_cost[(n,i)] = JuMP.@NLexpression(pm.model, dcline["cost"][1]*p_dc + dcline["cost"][2])
            elseif length(dcline["cost"]) == 3
                pmin = sum(dcline["pminf"][c] for c in _PM.conductor_ids(pm, n))
                pmax = sum(dcline["pmaxf"][c] for c in _PM.conductor_ids(pm, n))

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

                if report
                    sol(pm, n, :dcline, i)[:p_dc_sqr] = p_dc_sqr
                end

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
        for (n, nw_ref) in _PM.nws(pm))
    )
end


function _PM.constraint_model_voltage(pm::SOCWROAPowerModel, n::Int)
    w  = var(pm, n,  :w)
    wr = var(pm, n, :wr)
    wi = var(pm, n, :wi)


    for (i,j) in _PM.ids(pm, n, :buspairs)
        @NLconstraint(pm.model, (wr[(i,j)]^2 + wi[(i,j)]^2)/w[j] <= w[i])
    end
end

function _PM.constraint_thermal_limit_from(pm::SOCWROAPowerModel, n::Int, f_idx, rate_a)
    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    @NLconstraint(pm.model, sqrt(p_fr^2 + q_fr^2) <= rate_a)
end

function _PM.constraint_thermal_limit_to(pm::SOCWROAPowerModel, n::Int, t_idx, rate_a)
    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)
    @NLconstraint(pm.model, sqrt(p_to^2 + q_to^2) <= rate_a)
end



# Defines a variant of the SOCWRForm which forces the NL solver path

export NLSOCWRPowerModel

mutable struct NLSOCWRPowerModel <: _PM.AbstractSOCWRModel _PM.@pm_fields end



# Defines a variant of the QCWRTriForm without the linking constraints
export QCLSNoLinkPowerModel

mutable struct QCLSNoLinkPowerModel <: _PM.AbstractQCLSModel _PM.@pm_fields end


function _PM.constraint_model_voltage(pm::QCLSNoLinkPowerModel, n::Int)
    v = var(pm, n, :vm)
    t = var(pm, n, :va)

    td = var(pm, n, :td)
    si = var(pm, n, :si)
    cs = var(pm, n, :cs)

    w = var(pm, n, :w)
    wr = var(pm, n, :wr)
    lambda_wr = var(pm, n, :lambda_wr)
    wi = var(pm, n, :wi)
    lambda_wi = var(pm, n, :lambda_wi)

    for (i,b) in ref(pm, n, :bus)
        _IM.relaxation_sqr(pm.model, v[i], w[i])
    end

    for bp in _PM.ids(pm, n, :buspairs)
        i,j = bp
        @constraint(pm.model, t[i] - t[j] == td[bp])

        _PM.relaxation_sin(pm.model, td[bp], si[bp])
        _PM.relaxation_cos(pm.model, td[bp], cs[bp])
        _IM.relaxation_trilinear(pm.model, v[i], v[j], cs[bp], wr[bp], lambda_wr[bp,:])
        _IM.relaxation_trilinear(pm.model, v[i], v[j], si[bp], wi[bp], lambda_wi[bp,:])

        # this constraint is redudant and useful for debugging
        #_IM.relaxation_complex_product(pm.model, w[i], w[j], wr[bp], wi[bp])
   end

   for (i,branch) in ref(pm, n, :branch)
        pair = (branch["f_bus"], branch["t_bus"])
        buspair = ref(pm, n, :buspairs, pair)

        # to prevent this constraint from being posted on multiple parallel branchs
        if buspair["branch"] == i
            _PM.constraint_power_magnitude_sqr(pm, i, nw=n)
            _PM.constraint_power_magnitude_link(pm, i, nw=n)
        end
    end

end
