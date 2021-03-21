# MultiObjectiveProblems

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://manuelbb-upb.github.io/MultiObjectiveProblems.jl/dev)

This Julia package provides some test problems for Multiobjective Optimization.
It is an early work-in-progress side project.
Currently, only some real-valued problems with real-valued arguments are included,
mainly problems from the ZDT and DTLZ family.

Please refer to the [documentation](https://manuelbb-upb.github.io/MultiObjectiveProblems.jl/dev) for usage information.

## Installation
This package is not yet registered in the General Registry.
You can add it via 
```julia 
using Pkg
Pkg.add(url = "https://github.com/manuelbb-upb/MultiObjectiveProblems.jl.git")
```

## What about the heavy dependencies? `ForwardDiff`, `OSQP`, `NLopt` etc.

* `ForwardDiff` is used to provide gradient information for the objectives by means of automatic differentiation.
* `JuMP` and `OSQP` are required to calculate the *criticality* of a decision vector using the above gradients.
* `NLopt` is employed to check local optimality for multi-modal problems. As soon as Julia supports optional dependencies, `NLopt` will become one.
