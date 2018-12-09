function net2pmc_gen(net::Network; baseMVA=100, id_offset=0)
    gen_dict = Dict{String, Any}()
    to_warn = Int[]
    for r in eachrow(gen(net))
        old = Dict()
        if in("gen", keys(pmc(net))) && string(r[:element_id]) in keys(pmc(net)["gen"])
            old = pmc(net)["gen"][string(r[:element_id])]
        end
        there = !isempty(old)
        coeffs = there ? old["cost"] : [0.0]
        coeffs = !ismissing(r[:cost]) ? costcurve2pmc(cost_gen(net)[:coeffs][r[:cost]]) : coeffs
        gen_dict[string(r[:element_id] + id_offset)] = Dict(
            "index" => (!there || !ismissing(r[:element_id])) ? r[:element_id] : old["index"],
            "gen_bus" => (!there || !ismissing(r[:bus])) ? r[:bus] : old["gen_bus"],
            "pg" => (!there || !ismissing(r[:gen_p])) ? r[:gen_p] : old["pg"],
            "pmax" => (!there || !ismissing(r[:p_max])) ? r[:p_max] : old["pmax"],
            "pmin" => (!there || !ismissing(r[:p_min])) ? r[:p_min] : old["pmin"],
            "qg" => (!there || !ismissing(r[:gen_q])) ? r[:gen_q] : old["qg"],
            "qmax" => (!there || !ismissing(r[:q_max])) ? r[:q_max] : old["qmax"],
            "qmin" => (!there || !ismissing(r[:q_min])) ? r[:q_min] : old["qmin"],
            "startup" => (!there || !ismissing(r[:startup_cost])) ? r[:startup_cost] : old["startup"],
            "gen_status" => (!there || !ismissing(r[:status])) ? r[:status] : old["gen_status"],
            "ramp_10" => (!there || !ismissing(r[:ramp])) ? r[:ramp] : old["ramp_10"],
            "cost" => (!there || !ismissing(r[:cost]) ) ? coeffs : old["cost"], # Default or old value
            "model" => !there  ? -1 : old["model"], # Default or old value
            "ncost" => !there  ? 0 : old["ncost"], # Default or old value
            "qc1max" => 0.0,
            "qc2max" => 0.0,
            "mbase" => baseMVA,
            "pc2" => 0.0,
            "pc1" => 0.0,
            "ramp_q" => 0.0,
            "ramp_30" => 0.0,
            "apf" => 0.0,
            "shutdown" => 0.0,
            "ramp_agc" => 0.0,
            "vg" => 1.01,
            "qc1min" => 0.0,
            "qc2min" => 0.0,
        )
        if !there && isa(coeffs, PolynomialCost)
            gen_dict["cost"] = coefficients(coeffs)
            gen_dict["model"] = 2
            gen_dict["ncost"] = n_cost(coeffs)
        elseif !there && isa(coeffs, PWLCost)
            mw_pts = ustrip.(mws(coeffs))
            cost_pts = ustrip.(costs(coeffs))
            mix = vec(hcat(mw_pts, cost_pts)') # This is the format used by MATPOWER and PowerModels
            gen_dict["cost"] = mix
            gen_dict["model"] = 1
            gen_dict["ncost"] = length(mw_pts)
        else
            push!(to_warn, r[:element_id])
        end
    end
    if !isempty(to_warn)
        warn(LOGGER, "Using default costs for generator cost ids $to_warn")
    end
    return gen_dict
end

function net2pmc_ps_load(net::Network; baseMVA=100, id_offset=0)
    ps_load_dict = Dict{String, Any}()
    # Here we also need to add the price sensitive loads as generators
    to_warn = Int[]
    for r in eachrow(ps_load(net))
        old = Dict()
        if string(r[:element_id]) in keys(pmc(net)["gen"])
            old = pmc(net)["gen"][string(r[:element_id])]
        end
        there = !isempty(old)
        coeffs = there ? old["cost"] : [0.0]
        coeffs = !ismissing(r[:cost]) ? cost_load(net)[:coeffs][r[:cost]] : coeffs
        ps_load_dict[string(r[:element_id] + id_offset)] = Dict(
            "index" => (!there || !ismissing(r[:element_id])) ? r[:element_id] + id_offset : old["index"],
            "gen_bus" => (!there || !ismissing(r[:bus])) ? r[:bus] : old["gen_bus"],
            "pmax" => (!there || !ismissing(r[:load_max])) ? r[:load_max] : old["pmax"],
            "cost" => (!there || !ismissing(r[:cost]) ) ? costcurve2pmc(cost_load(net)[:coeffs][r[:cost]]) : old["cost"], # Default or old value
            "model" => !there  ? -1 : old["model"], # Default or old value
            "ncost" => !there  ? 0 : old["ncost"], # Default or old value
            "gen_status" => (!there || !ismissing(r[:status])) ? r[:status] : old["gen_status"],
            "ramp_10" => 0.0,
            "startup" => 0.0,
            "qg" => 0.0,
            "qmax" => 0.0,
            "qmin" => 0.0,
            "pmin" => 0.0,
            "qc1max" => 0.0,
            "model" => 2,
            "qc2max" => 0.0,
            "mbase" => baseMVA,
            "pc2" => 0.0,
            "pc1" => 0.0,
            "ramp_q" => 0.0,
            "ramp_30" => 0.0,
            "apf" => 0.0,
            "shutdown" => 0.0,
            "ramp_agc" => 0.0,
            "vg" => 1.01,
            "qc1min" => 0.0,
            "qc2min" => 0.0,
        )
        # TODO: Check the signs of ps load costs
        if !there && isa(coeffs, PolynomialCost)
            # Load is negative generation. By serving load, the total cost receives a negative contribution
            # due to the amount payed by the utility.
            ps_load_dict["cost"] = -coefficients(coeffs) # The cost curve has to be flipped along the x axis.
            aux = ps_load_dict["pmin"]
            ps_load_dict["pmin"] = -gen_dict["pmax"] # The cost curve has to be flipped along the y axis.
            ps_load_dict["pmax"] = -aux
            ps_load_dict["model"] = 2
            ps_load_dict["ncost"] = n_cost(coeffs)
        elseif !there && isa(coeffs, PWLCost)
            mw_pts = -ustrip.(mws(coeffs)) # The mws have to be flipped.
            cost_pts = -ustrip.(costs(coeffs)) # The cost has to be flipped too. It's a cost for the utility, profit for generators.
            mix = vec(hcat(mw_pts, cost_pts)') # This is the format used by MATPOWER and PowerModels
            ps_load_dict["cost"] = mix
            ps_load_dict["model"] = 1
            ps_load_dict["ncost"] = length(mw_pts)
            ps_load_dict["pmin"] = -ps_load_dict["pmax"] # The cost curve has to be flipped along the y axis. Load is negative generation
            ps_load_dict["pmax"] = 0.0
        else
            push!(to_warn, r[:element_id])
        end
    end
    if !isempty(to_warn)
        warn(LOGGER, "Used default costs for ps load cost ids $to_warn.")
    end
    return ps_load_dict
end
