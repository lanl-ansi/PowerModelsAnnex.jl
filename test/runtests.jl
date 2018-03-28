using PowerModelsAnnex
using Unitful
using PowerModelsAnnex.Units
PMA = PowerModelsAnnex

using Logging
# suppress warnings during testing
Logging.configure(level=ERROR)

using Memento
setlevel!(getlogger(InfrastructureModels), "error")
setlevel!(getlogger(PowerModels), "error")

using JuMP
using PowerModels
PMs = PowerModels

using Ipopt

using Base.Test

# default setup for solvers
ipopt_solver = IpoptSolver(tol=1e-6, print_level=0)

# this will work because PowerModels is a dependency
pmpath = Pkg.dir("PowerModels") # doing like this so that it does not depend on relative paths
casepath = joinpath(pmpath, "test", "data")
case_files = [
    "../../PowerModels/test/data/matpower/case3.m",
    "../../PowerModels/test/data/matpower/case5.m",
    "../../PowerModels/test/data/matpower/case5_asym.m",
    "../../PowerModels/test/data/matpower/case5_dc.m",
    "../../PowerModels/test/data/matpower/case14.m",
    "../../PowerModels/test/data/matpower/case24.m",
    "../../PowerModels/test/data/matpower/case30.m"
]

@testset "PowerModelsAnnex" begin

include("form/wr.jl")

include("pglib/api.jl")
include("pglib/sad.jl")

include("model/pf.jl")
include("model/opf.jl")

include("frontend/frontend.jl")

end
