PowerModelsAnnex.jl Change Log
==============================

### Staged
- nothing

### v0.11.0
- Update PowerModels v0.21 and JuMP's new nonlinear interface (breaking)

### v0.10.0
- Update PowerModels v0.20
- Drop `frontend` directory, lack of maintenance on dependent packages
- Drop `islanding` directory, no tests

### v0.9.0
- Revise `solve_opf_api` to consider reactive power dispatch

### v0.8.5
- Add closest operating point formulation `solve_opf_cop`

### v0.8.4
- Update PowerModels `run_*` functions to `solve_*`
- Add support for Memento v1.4

### v0.8.3
- Add support for JuMP v1.0

### v0.8.2
- Add support for JuMP v0.23
- Update minimum Julia version to v1.6 (LTS)

### v0.8.1
- Add support for Memento v1.3

### v0.8.0
- Update to JuMP v0.22, PowerModels v0.19
- Drop support for JuMP v0.21
- Remove dependency on MathOptInterface package

### v0.7.1
- Add support for Memento v1.2
- Update use of `with_optimizer` to `optimizer_with_attributes`

### v0.7.0
- Update to InfrastructureModels v0.6 and PowerModels v0.18

### v0.6.1
- Add variants of piecewise linear cost model formulations
- Update to DataFrames v0.21

### v0.6.0
- Update to PowerModels v0.17
- Added support for Memento v1.1

### v0.5.0
- Update to PowerModels v0.16

### v0.4.4
- Fixed solution reporting bug in api model

### v0.4.3
- Update pacakge internal short names to `_PM` and `_IM`
- Add support for optional branch power flows in api and sad models

### v0.4.2
- Add support for Memento v0.13, v1.0

### v0.4.1
- Add support for JuMP v0.21
- Drop Manifest.toml (#57)

### v0.4.0
- Update to PowerModels v0.15

### v0.3.1
- Update to InfrastructureModels v0.4

### v0.3.0
- Update to PowerModels v0.14 and new naming conventions (#65)

### v0.2.6
- Updates for PowerModels v0.13

### v0.2.5
- Updates to frontend for MOI status values (#60)

### v0.2.4
- Updates for JuMP v0.20 and Julia v1.2 (#59)
- Resolve Dataframes deprecations (#58)

### v0.2.3
- Updates for PowerSystemsUnits
- Fixes for Frontend (#53)

### v0.2.2
- Update to PowerModels v0.12

### v0.2.1
- Update to PowerModels v0.11

### v0.2.0
- Update to JuMP v0.19/MathOptInterface

### v0.1.12
- Dropped support for Julia v0.6/v0.7

### Previous
- See Github releases for details
