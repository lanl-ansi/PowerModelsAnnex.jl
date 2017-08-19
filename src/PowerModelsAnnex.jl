isdefined(Base, :__precompile__) && __precompile__()

module PowerModelsAnnex

using Compat
using JuMP
using PowerModels
PMs = PowerModels

include("form/wr.jl")

include("model/pf.jl")
include("model/opf.jl")

end
