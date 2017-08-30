isdefined(Base, :__precompile__) && __precompile__()

module PowerModelsAnnex

using Compat
using JuMP
using PowerModels
PMs = PowerModels

include("form/wr.jl")

include("model/pf.jl")
include("model/opf.jl")

include("pglib/shared.jl")
include("pglib/api.jl")
include("pglib/sad.jl")

end
