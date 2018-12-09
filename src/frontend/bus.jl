function net2pmc_bus(net::Network)
    bus_dict = Dict{String, Any}()
    for r in eachrow(bus(net))
        old = Dict()
        if in("bus", keys(pmc(net))) && string(r[:element_id]) in keys(pmc(net)["bus"])
            old = pmc(net)["bus"][string(r[:element_id])]
        end
        there = !isempty(old)
        bus_dict[string(r[:element_id])] = Dict(
            "index" => (!there || !ismissing(r[:element_id])) ? r[:element_id] : old["index"],
            "bus_name" => (!there || !ismissing(r[:element_id])) ? r[:name] : old["bus_name"],
            "vm" => (!there || !ismissing(r[:element_id])) ? r[:voltage] : old["vm"],
            "vmax" => (!there || !ismissing(r[:element_id])) ? r[:volt_max] : old["vmax"],
            "vmin" => (!there || !ismissing(r[:element_id])) ? r[:volt_min] : old["vmin"],
            "bus_type" => (!there || !ismissing(r[:element_id])) ? r[:bus_type] : old["bus_type"],
            "base_kv" => (!there || !ismissing(r[:element_id])) ? r[:base_kv] : old["base_kv"],
            "zone" => (!there || !ismissing(r[:element_id])) ? r[:zone] : old["zone"],
            "bus_i" => (!there || !ismissing(r[:element_id])) ? r[:element_id] : old["bus_i"],
            "area" => (!there || !ismissing(r[:element_id])) ? r[:area] : old["area"],
            "va" => 0.0,
        )
    end
    return bus_dict
end
