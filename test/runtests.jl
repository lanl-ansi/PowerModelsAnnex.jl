using PowerModelsAnnex

using Logging
# suppress warnings during testing
Logging.configure(level=ERROR)

using JuMP
using PowerModels
PMs = PowerModels

using Ipopt

using Base.Test

# default setup for solvers
ipopt_solver = IpoptSolver(tol=1e-6, print_level=0)


# this will work because PowerModels is a dependency
case_files = [
    "../../PowerModels/test/data/case2.m",
    "../../PowerModels/test/data/case3.m",
    "../../PowerModels/test/data/case5.m",
    "../../PowerModels/test/data/case5_asym.m",
    "../../PowerModels/test/data/case5_dc.m",
    "../../PowerModels/test/data/case14.m",
    "../../PowerModels/test/data/case24.m",
    "../../PowerModels/test/data/case30.m"
]

include("model/pf.jl")
include("model/opf.jl")

