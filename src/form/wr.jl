

# Defines a variant of the SOCWRForm which is appropriate for outer approximation

export SOCWROAPowerModel

mutable struct SOCWROAPowerModel <: _PM.AbstractSOCWRModel _PM.@pm_fields end


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
