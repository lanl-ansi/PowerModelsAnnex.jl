using PowerModelsAnnex

using Logging
# suppress warnings during testing
Logging.configure(level=ERROR)

using Ipopt

using Base.Test

# default setup for solvers
ipopt_solver = IpoptSolver(tol=1e-6, print_level=0)

