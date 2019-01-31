export post_ac_opf, post_soc_opf, post_qc_opf, post_dc_opf

"""
Given a JuMP model and a PowerModels network data structure,
Builds an AC-OPF formulation of the given data and returns the JuMP model
"""
function post_ac_opf(data::Dict{String,Any}, model=JuMP.Model())
    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    PowerModels.standardize_cost_terms(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]

    JuMP.@variable(model, va[i in keys(ref[:bus])])
    JuMP.@variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)

    JuMP.@variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    JuMP.@variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

    JuMP.@variable(model, ref[:arcs_dc_param][a]["pmin"] <= p_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["pmax"])
    JuMP.@variable(model, ref[:arcs_dc_param][a]["qmin"] <= q_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["qmax"])

    from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])
    JuMP.@objective(model, Min,
        sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) +
        sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
    )

    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        JuMP.@constraint(model, va[i] == 0)
    end

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Bus KCL
        JuMP.@constraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) +
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts)*vm[i]^2
        )
        JuMP.@constraint(model,
            sum(q[a] for a in ref[:bus_arcs][i]) +
            sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) -
            sum(load["qd"] for load in bus_loads) +
            sum(shunt["bs"] for shunt in bus_shunts)*vm[i]^2
        )
    end

    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        vm_fr = vm[branch["f_bus"]]
        vm_to = vm[branch["t_bus"]]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        # Line Flow
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2

        # AC Line Flow Constraints
        JuMP.@NLconstraint(model, p_fr ==  (g+g_fr)/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        JuMP.@NLconstraint(model, q_fr == -(b+b_fr)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        JuMP.@NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        JuMP.@NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )

        # Phase Angle Difference Limit
        JuMP.@constraint(model, va_fr - va_to <= branch["angmax"])
        JuMP.@constraint(model, va_fr - va_to >= branch["angmin"])

        # Apparent Power Limit, From and To
        JuMP.@constraint(model, p[f_idx]^2 + q[f_idx]^2 <= branch["rate_a"]^2)
        JuMP.@constraint(model, p[t_idx]^2 + q[t_idx]^2 <= branch["rate_a"]^2)
    end

    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])

        JuMP.@constraint(model, (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    end

    return model
end


"""
Given a JuMP model and a PowerModels network data structure,
Builds an SOC-OPF formulation of the given data and returns the JuMP model
"""
function post_soc_opf(data::Dict{String,Any}, model=JuMP.Model())
    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    PowerModels.standardize_cost_terms(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]

    JuMP.@variable(model, ref[:bus][i]["vmin"]^2 <= w[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"]^2, start=1.001)

    wr_min, wr_max, wi_min, wi_max = PowerModels.calc_voltage_product_bounds(ref[:buspairs])

    JuMP.@variable(model, wr_min[bp] <= wr[bp in keys(ref[:buspairs])] <= wr_max[bp], start=1.0)
    JuMP.@variable(model, wi_min[bp] <= wi[bp in keys(ref[:buspairs])] <= wi_max[bp])

    JuMP.@variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    JuMP.@variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

    JuMP.@variable(model, ref[:arcs_dc_param][a]["pmin"] <= p_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["pmax"])
    JuMP.@variable(model, ref[:arcs_dc_param][a]["qmin"] <= q_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["qmax"])

    from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])
    JuMP.@objective(model, Min,
        sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) +
        sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
    )

    for (bp, buspair) in ref[:buspairs]
        i,j = bp

        # Voltage Product Relaxation Lowerbound
        JuMP.@constraint(model, wr[(i,j)]^2 + wi[(i,j)]^2 <= w[i]*w[j])

        vfub = buspair["vm_fr_max"]
        vflb = buspair["vm_fr_min"]
        vtub = buspair["vm_to_max"]
        vtlb = buspair["vm_to_min"]
        tdub = buspair["angmax"]
        tdlb = buspair["angmin"]

        phi = (tdub + tdlb)/2
        d   = (tdub - tdlb)/2

        sf = vflb + vfub
        st = vtlb + vtub

        # Voltage Product Relaxation Upperbound
        JuMP.@constraint(model, sf*st*(cos(phi)*wr[(i,j)] + sin(phi)*wi[(i,j)]) - vtub*cos(d)*st*w[i] - vfub*cos(d)*sf*w[j] >=  vfub*vtub*cos(d)*(vflb*vtlb - vfub*vtub))
        JuMP.@constraint(model, sf*st*(cos(phi)*wr[(i,j)] + sin(phi)*wi[(i,j)]) - vtlb*cos(d)*st*w[i] - vflb*cos(d)*sf*w[j] >= -vflb*vtlb*cos(d)*(vflb*vtlb - vfub*vtub))
    end

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Bus KCL
        JuMP.@constraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) +
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts)*w[i]
        )
        JuMP.@constraint(model,
            sum(q[a] for a in ref[:bus_arcs][i]) +
            sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) -
            sum(load["qd"] for load in bus_loads) +
            sum(shunt["bs"] for shunt in bus_shunts)*w[i]
        )
    end

    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])
        bp_idx = (branch["f_bus"], branch["t_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        w_fr = w[branch["f_bus"]]
        w_to = w[branch["t_bus"]]
        wr_br = wr[bp_idx]
        wi_br = wi[bp_idx]

        # Line Flow
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2

        # AC Line Flow Constraints
        JuMP.@constraint(model, p_fr ==  (g+g_fr)/tm*w_fr + (-g*tr+b*ti)/tm*(wr_br) + (-b*tr-g*ti)/tm*(wi_br) )
        JuMP.@constraint(model, q_fr == -(b+b_fr)/tm*w_fr - (-b*tr-g*ti)/tm*(wr_br) + (-g*tr+b*ti)/tm*(wi_br) )

        JuMP.@constraint(model, p_to ==  (g+g_to)*w_to + (-g*tr-b*ti)/tm*(wr_br) + (-b*tr+g*ti)/tm*(-wi_br) )
        JuMP.@constraint(model, q_to == -(b+b_to)*w_to - (-b*tr+g*ti)/tm*(wr_br) + (-g*tr-b*ti)/tm*(-wi_br) )

        # Phase Angle Difference Limit
        JuMP.@constraint(model, wi_br <= tan(branch["angmax"])*wr_br)
        JuMP.@constraint(model, wi_br >= tan(branch["angmin"])*wr_br)

        # Apparent Power Limit, From and To
        JuMP.@constraint(model, p[f_idx]^2 + q[f_idx]^2 <= branch["rate_a"]^2)
        JuMP.@constraint(model, p[t_idx]^2 + q[t_idx]^2 <= branch["rate_a"]^2)
    end

    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])

        JuMP.@constraint(model, (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    end

    return model
end


"""
Given a JuMP model and a PowerModels network data structure,
Builds an QC-OPF formulation of the given data and returns the JuMP model
Implementation provided by @sidhant172
"""
function post_qc_opf(data::Dict{String,Any}, model=JuMP.Model())
    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    PowerModels.standardize_cost_terms(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]

    # voltage angle and magnitude
    JuMP.@variable(model, va[i in keys(ref[:bus])])
    JuMP.@variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)

    # voltage squared
    JuMP.@variable(model, ref[:bus][i]["vmin"]^2 <= w[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"]^2, start=1.001)

    # voltage product
    wr_min, wr_max, wi_min, wi_max = PowerModels.calc_voltage_product_bounds(ref[:buspairs])
    JuMP.@variable(model, wr_min[bp] <= wr[bp in keys(ref[:buspairs])] <= wr_max[bp], start=1.0)
    JuMP.@variable(model, wi_min[bp] <= wi[bp in keys(ref[:buspairs])] <= wi_max[bp])

    # voltage angle differences
    JuMP.@variable(model, ref[:buspairs][bp]["angmin"] <= td[bp in keys(ref[:buspairs])]  <= ref[:buspairs][bp]["angmax"])

    # variable multipliers in lambda formulation
    JuMP.@variable(model, 0 <= lambda_wr[bp in keys(ref[:buspairs]), 1:8] <= 1)
    JuMP.@variable(model, 0 <= lambda_wi[bp in keys(ref[:buspairs]), 1:8] <= 1)


    # cosine variables

    # computing bounds for cosine variables
    cos_min = Dict([(bp, -Inf) for bp in keys(ref[:buspairs])])
    cos_max = Dict([(bp, Inf) for bp in keys(ref[:buspairs])])

    for (bp, buspair) in ref[:buspairs]
        angmin = buspair["angmin"]
        angmax = buspair["angmax"]
        if angmin >= 0
            cos_max[bp] = cos(angmin)
            cos_min[bp] = cos(angmax)
        end
        if angmax <= 0
            cos_max[bp] = cos(angmax)
            cos_min[bp] = cos(angmin)
        end
        if angmin < 0 && angmax > 0
            cos_max[bp] = 1.0
            cos_min[bp] = min(cos(angmin), cos(angmax))
        end
    end
    # end computing bounds for cosine variables

    # defining cosine variables
    JuMP.@variable(model,  cos_min[bp] <= cs[bp in keys(ref[:buspairs])] <= cos_max[bp])

    # defining sine variables
    JuMP.@variable(model, sin(ref[:buspairs][bp]["angmin"]) <= si[bp in keys(ref[:buspairs])] <=  sin(ref[:buspairs][bp]["angmax"]))

    # current magnitude squared
    # compute upper bound
    cm_ub = Dict()
    for (bp, buspair) in ref[:buspairs]
        cm_ub[bp] = ((buspair["rate_a"]*buspair["tap"])/buspair["vm_fr_min"])^2
    end
    # define current magnitude variable
    JuMP.@variable(model, 0 <= cm[bp in keys(ref[:buspairs])] <= cm_ub[bp])


    # line flow variables
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

    # generation pg an qg
    JuMP.@variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    JuMP.@variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    # dc line flows
    JuMP.@variable(model, ref[:arcs_dc_param][a]["pmin"] <= p_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["pmax"])
    JuMP.@variable(model, ref[:arcs_dc_param][a]["qmin"] <= q_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["qmax"])


    # objective
    from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])
    JuMP.@objective(model, Min,
        sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) +
        sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
    )



    #### beginning of constraints ####

    # relaxation of vm square
    for (i,bus) in ref[:bus]
        JuMP.@constraint(model, w[i] >= vm[i]^2)
        JuMP.@constraint(model, w[i] <= (bus["vmin"] + bus["vmax"])*vm[i] - bus["vmin"]*bus["vmax"])
    end

    # buspair voltage constraints
    for (bp, buspair) in ref[:buspairs]
        i,j = bp
        JuMP.@constraint(model, va[i] - va[j] == td[bp])

        # relaxation sin
        ub = buspair["angmax"]
        lb = buspair["angmin"]
        max_ad = max(abs(lb),abs(ub))
        if lb < 0 && ub > 0
            JuMP.@constraint(model, si[bp] <= cos(max_ad/2)*(td[bp] - max_ad/2) + sin(max_ad/2))
            JuMP.@constraint(model, si[bp] >= cos(max_ad/2)*(td[bp] + max_ad/2) - sin(max_ad/2))
        end
        if ub <= 0
            JuMP.@constraint(model, si[bp] <= (sin(lb) - sin(ub))/(lb-ub)*(td[bp] - lb) + sin(lb))
            JuMP.@constraint(model, si[bp] >= cos(max_ad/2)*(td[bp] + max_ad/2) - sin(max_ad/2))
        end
        if lb >= 0
            JuMP.@constraint(model, si[bp] <= cos(max_ad/2)*(td[bp] - max_ad/2) + sin(max_ad/2))
            JuMP.@constraint(model, si[bp] >= (sin(lb) - sin(ub))/(lb-ub)*(td[bp] - lb) + sin(lb))
        end
        # end of relaxation sin

        # relaxation cos
        JuMP.@constraint(model, cs[bp] <= 1 - (1-cos(max_ad))/(max_ad*max_ad)*(td[bp]^2))
        JuMP.@constraint(model, cs[bp] >= (cos(lb) - cos(ub))/(lb-ub)*(td[bp] - lb) + cos(lb))
        # end of relaxation cos

        ##### relaxation trilinear wr
        wr_val = [
            ref[:bus][i]["vmin"] * ref[:bus][j]["vmin"] * cos_min[bp]
            ref[:bus][i]["vmin"] * ref[:bus][j]["vmin"] * cos_max[bp]
            ref[:bus][i]["vmin"] * ref[:bus][j]["vmax"] * cos_min[bp]
            ref[:bus][i]["vmin"] * ref[:bus][j]["vmax"] * cos_max[bp]
            ref[:bus][i]["vmax"] * ref[:bus][j]["vmin"] * cos_min[bp]
            ref[:bus][i]["vmax"] * ref[:bus][j]["vmin"] * cos_max[bp]
            ref[:bus][i]["vmax"] * ref[:bus][j]["vmax"] * cos_min[bp]
            ref[:bus][i]["vmax"] * ref[:bus][j]["vmax"] * cos_max[bp]
        ]

        JuMP.@constraint(model, wr[bp] == sum(wr_val[ii]*lambda_wr[bp,ii] for ii in 1:8))
        JuMP.@constraint(model, vm[i] == (lambda_wr[bp,1] + lambda_wr[bp,2] + lambda_wr[bp,3] + lambda_wr[bp,4])*ref[:bus][i]["vmin"] +
                            (lambda_wr[bp,5] + lambda_wr[bp,6] + lambda_wr[bp,7] + lambda_wr[bp,8])*ref[:bus][i]["vmax"])
        JuMP.@constraint(model, vm[j] == (lambda_wr[bp,1] + lambda_wr[bp,2] + lambda_wr[bp,5] + lambda_wr[bp,6])*ref[:bus][j]["vmin"] +
                            (lambda_wr[bp,3] + lambda_wr[bp,4] + lambda_wr[bp,7] + lambda_wr[bp,8])*ref[:bus][j]["vmax"])
        JuMP.@constraint(model, cs[bp] == (lambda_wr[bp,1] + lambda_wr[bp,3] + lambda_wr[bp,5] + lambda_wr[bp,7])*cos_min[bp] +
                            (lambda_wr[bp,2] + lambda_wr[bp,4] + lambda_wr[bp,6] + lambda_wr[bp,8])*cos_max[bp])
        JuMP.@constraint(model, sum(lambda_wr[bp,ii] for ii in 1:8) == 1)
        #### end of relaxation trilinear wr

        ##### relaxation trilinear wi
        wi_val = [
            ref[:bus][i]["vmin"] * ref[:bus][j]["vmin"] * sin(buspair["angmin"])
            ref[:bus][i]["vmin"] * ref[:bus][j]["vmin"] * sin(buspair["angmax"])
            ref[:bus][i]["vmin"] * ref[:bus][j]["vmax"] * sin(buspair["angmin"])
            ref[:bus][i]["vmin"] * ref[:bus][j]["vmax"] * sin(buspair["angmax"])
            ref[:bus][i]["vmax"] * ref[:bus][j]["vmin"] * sin(buspair["angmin"])
            ref[:bus][i]["vmax"] * ref[:bus][j]["vmin"] * sin(buspair["angmax"])
            ref[:bus][i]["vmax"] * ref[:bus][j]["vmax"] * sin(buspair["angmin"])
            ref[:bus][i]["vmax"] * ref[:bus][j]["vmax"] * sin(buspair["angmax"])
        ]

        JuMP.@constraint(model, wi[bp] == sum(wi_val[ii]*lambda_wi[bp,ii] for ii in 1:8))
        JuMP.@constraint(model, vm[i] == (lambda_wi[bp,1] + lambda_wi[bp,2] + lambda_wi[bp,3] + lambda_wi[bp,4])*ref[:bus][i]["vmin"] +
                            (lambda_wi[bp,5] + lambda_wi[bp,6] + lambda_wi[bp,7] + lambda_wi[bp,8])*ref[:bus][i]["vmax"])
        JuMP.@constraint(model, vm[j] == (lambda_wi[bp,1] + lambda_wi[bp,2] + lambda_wi[bp,5] + lambda_wi[bp,6])*ref[:bus][j]["vmin"] +
                            (lambda_wi[bp,3] + lambda_wi[bp,4] + lambda_wi[bp,7] + lambda_wi[bp,8])*ref[:bus][j]["vmax"])
        JuMP.@constraint(model, si[bp] == (lambda_wi[bp,1] + lambda_wi[bp,3] + lambda_wi[bp,5] + lambda_wi[bp,7])*sin(buspair["angmin"]) +
                            (lambda_wi[bp,2] + lambda_wi[bp,4] + lambda_wi[bp,6] + lambda_wi[bp,8])*sin(buspair["angmax"]))
        JuMP.@constraint(model, sum(lambda_wi[bp,ii] for ii in 1:8) == 1)
        #### end of relaxation trilinear wi

        # relaxation tighten vv
        val = [
            ref[:bus][i]["vmin"] * ref[:bus][j]["vmin"]
            ref[:bus][i]["vmin"] * ref[:bus][j]["vmin"]
            ref[:bus][i]["vmin"] * ref[:bus][j]["vmax"]
            ref[:bus][i]["vmin"] * ref[:bus][j]["vmax"]
            ref[:bus][i]["vmax"] * ref[:bus][j]["vmin"]
            ref[:bus][i]["vmax"] * ref[:bus][j]["vmin"]
            ref[:bus][i]["vmax"] * ref[:bus][j]["vmax"]
            ref[:bus][i]["vmax"] * ref[:bus][j]["vmax"]
        ]

        JuMP.@constraint(model, sum(lambda_wr[bp,ii]*val[ii] - lambda_wi[bp,ii]*val[ii] for ii in 1:8) == 0)
        # end of relaxation tighten vv

        vfub = buspair["vm_fr_max"]
        vflb = buspair["vm_fr_min"]
        vtub = buspair["vm_to_max"]
        vtlb = buspair["vm_to_min"]
        tdub = buspair["angmax"]
        tdlb = buspair["angmin"]

        phi = (tdub + tdlb)/2
        d   = (tdub - tdlb)/2

        sf = vflb + vfub
        st = vtlb + vtub

        # Voltage Product Relaxation Upperbound
        JuMP.@constraint(model, sf*st*(cos(phi)*wr[(i,j)] + sin(phi)*wi[(i,j)]) - vtub*cos(d)*st*w[i] - vfub*cos(d)*sf*w[j] >=  vfub*vtub*cos(d)*(vflb*vtlb - vfub*vtub))
        JuMP.@constraint(model, sf*st*(cos(phi)*wr[(i,j)] + sin(phi)*wi[(i,j)]) - vtlb*cos(d)*st*w[i] - vflb*cos(d)*sf*w[j] >= -vflb*vtlb*cos(d)*(vflb*vtlb - vfub*vtub))
    end


   for (i,branch) in ref[:branch]
        bp = (branch["f_bus"], branch["t_bus"])
        buspair = ref[:buspairs][bp]
        tm = branch["tap"]

        # to prevent this constraint from being posted on multiple parallel branches
        if buspair["branch"] == i
            # extract quantities
            g, b = PowerModels.calc_branch_y(branch)
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"]
            b_fr = branch["b_fr"]


            # extract variables
            p_fr = p[(i,branch["f_bus"],branch["t_bus"])]
            q_fr = q[(i,branch["f_bus"],branch["t_bus"])]
            w_fr = w[branch["f_bus"]]
            w_to = w[branch["t_bus"]]

            JuMP.@constraint(model, p_fr^2 + q_fr^2 <= w_fr/tm^2*cm[bp])

            ym_sh_sqr = g_fr^2 + b_fr^2

            JuMP.@constraint(model, cm[bp] == (g^2 + b^2)*(w_fr/tm^2 + w_to - 2*(tr*wr[bp] + ti*wi[bp])/tm^2) - ym_sh_sqr*(w_fr/tm^2) + 2*(g_fr*p_fr - b_fr*q_fr))
        end
    end

    # end of QC tri-form voltage constraint
    # end of constraint voltage


    # constraint theta ref
    for i in keys(ref[:ref_buses])
        JuMP.@constraint(model, va[i] == 0)
    end

    # constraint KCL shunt
    for (i,bus) in ref[:bus]
        # Bus KCL
        gs = 0.0
        bs = 0.0
        if length(ref[:bus_shunts][i]) > 0
            shunt_num = ref[:bus_shunts][i][1]
            gs = ref[:shunt][shunt_num]["gs"]
            bs = ref[:shunt][shunt_num]["bs"]
        end

        JuMP.@constraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) +
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) - sum(load["pd"] for (l,load) in ref[:load] if load["load_bus"]==i)  - gs*w[i]
        )
        JuMP.@constraint(model,
            sum(q[a] for a in ref[:bus_arcs][i]) +
            sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) - sum(load["qd"] for (l,load) in ref[:load] if load["load_bus"]==i) + bs*w[i]
        )
    end


    # ohms laws, voltage angle difference limits
    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])
        bp_idx = (branch["f_bus"], branch["t_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        w_fr = w[branch["f_bus"]]
        w_to = w[branch["t_bus"]]
        wr_br = wr[bp_idx]
        wi_br = wi[bp_idx]

        # Line Flow
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        g_to = branch["g_to"]
        b_fr = branch["b_fr"]
        b_to = branch["b_to"]
        tm = branch["tap"]

        # AC Line Flow Constraints
        JuMP.@constraint(model, p_fr ==  (g+g_fr)/tm^2*w_fr + (-g*tr+b*ti)/tm^2*wr_br + (-b*tr-g*ti)/tm^2*wi_br )
        JuMP.@constraint(model, q_fr == -(b+b_fr)/tm^2*w_fr - (-b*tr-g*ti)/tm^2*wr_br + (-g*tr+b*ti)/tm^2*wi_br )

        JuMP.@constraint(model, p_to ==  (g+g_to)*w_to + (-g*tr-b*ti)/tm^2*wr_br + (-b*tr+g*ti)/tm^2*-wi_br )
        JuMP.@constraint(model, q_to == -(b+b_to)*w_to - (-b*tr+g*ti)/tm^2*wr_br + (-g*tr-b*ti)/tm^2*-wi_br )

        # Phase Angle Difference Limit
        JuMP.@constraint(model, wi_br <= tan(branch["angmax"])*wr_br)
        JuMP.@constraint(model, wi_br >= tan(branch["angmin"])*wr_br)

        # Apparent Power Limit, From and To
        JuMP.@constraint(model, p[f_idx]^2 + q[f_idx]^2 <= branch["rate_a"]^2)
        JuMP.@constraint(model, p[t_idx]^2 + q[t_idx]^2 <= branch["rate_a"]^2)
    end

    # DC line constraints
    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])

        JuMP.@constraint(model, (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    end

    return model
end



"""
Given a JuMP model and a PowerModels network data structure,
Builds an DC-OPF formulation of the given data and returns the JuMP model
"""
function post_dc_opf(data::Dict{String,Any}, model=JuMP.Model())
    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    PowerModels.standardize_cost_terms(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]

    JuMP.@variable(model, va[i in keys(ref[:bus])])

    JuMP.@variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])

    JuMP.@variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs_from]] <= ref[:branch][l]["rate_a"])

    p_expr = Dict([((l,i,j), 1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]])
    p_expr = merge(p_expr, Dict([((l,j,i), -1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]]))

    JuMP.@variable(model, ref[:arcs_dc_param][a]["pmin"] <= p_dc[a in ref[:arcs_dc]] <= ref[:arcs_dc_param][a]["pmax"])

    from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])
    JuMP.@objective(model, Min,
        sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) +
        sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
    )

    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        JuMP.@constraint(model, va[i] == 0)
    end

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Bus KCL
        JuMP.@constraint(model,
            sum(p_expr[a] for a in ref[:bus_arcs][i]) +
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts)*1.0^2
        )
    end

    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])

        p_fr = p[f_idx]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        g, b = PowerModels.calc_branch_y(branch)

        # DC Line Flow Constraints
        JuMP.@constraint(model, p_fr == -b*(va_fr - va_to))

        # Phase Angle Difference Limit
        JuMP.@constraint(model, va_fr - va_to <= branch["angmax"])
        JuMP.@constraint(model, va_fr - va_to >= branch["angmin"])

        # Apparent Power Limit, From and To
        # covered by variable bounds
    end

    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])

        JuMP.@constraint(model, (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    end

    return model
end

