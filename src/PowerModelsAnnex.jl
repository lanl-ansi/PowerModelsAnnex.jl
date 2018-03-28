isdefined(Base, :__precompile__) && __precompile__()

module PowerModelsAnnex

using Compat
using JuMP
using PowerModels
PMs = PowerModels

include("form/acr.jl")
include("form/wr.jl")
include("form/shared.jl")

include("model/pf.jl")
include("model/opf.jl")

include("pglib/shared.jl")
include("pglib/api.jl")
include("pglib/sad.jl")

include("frontend/frontend.jl")
using PowerModelsAnnex.FrontEnd

end
