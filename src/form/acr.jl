export NLACRPowerModel, NLACRForm

abstract type NLACRForm <: PMs.AbstractACRForm end

const NLACRPowerModel = PMs.GenericPowerModel{NLACRForm}

"default NLACRForm constructor"
NLACRPowerModel(data::Dict{String,Any}; kwargs...) = PMs.GenericPowerModel(data, NLACRForm; kwargs...)