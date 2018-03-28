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
pmpath = Pkg.dir("PowerModels") # doing like this so that it does not depend on relative paths
casepath = joinpath(pmpath, "test", "data")
case_files = [
    joinpath(casepath, "case3.m"),
    joinpath(casepath, "case5.m"),
    joinpath(casepath, "case5_asym.m"),
    joinpath(casepath, "case5_dc.m"),
    joinpath(casepath, "case14.m"),
    joinpath(casepath, "case24.m"),
    joinpath(casepath, "case30.m")
]

include("form/wr.jl")

include("pglib/api.jl")
include("pglib/sad.jl")

include("model/pf.jl")
include("model/opf.jl")

# include("frontend/frontend.jl")
