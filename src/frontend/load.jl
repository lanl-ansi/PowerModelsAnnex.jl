function net2pmc_load(net::Network)
    load_dict = Dict{String, Any}()
    for r in eachrow(pi_load(net))
        old = Dict()
        if in("load", keys(pmc(net))) && string(r[:element_id]) in keys(pmc(net)["load"])
            old = pmc(net)["load"][string(r[:element_id])]
        end
        there = !isempty(old)
        load_dict[string(r[:element_id])] = Dict(
            "index" => (!there || !ismissing(r[:element_id])) ? r[:element_id] : old["index"],
            "load_bus" => (!there || !ismissing(r[:element_id])) ? r[:bus] : old["load_bus"],
            "qd" => (!there || !ismissing(r[:element_id])) ? r[:load_q] : old["qd"],
            "pd" => (!there || !ismissing(r[:element_id])) ? r[:load_p] : old["pd"],
            "status" => (!there || !ismissing(r[:element_id])) ? r[:status] : old["status"],
        )
    end
    return load_dict
end
