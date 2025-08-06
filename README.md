# cuOpt.jl

[cuOpt.jl](https://github.com/jump-dev/cuOpt.jl) is a wrapper for [NVIDIA cuOpt](https://github.com/NVIDIA/cuOpt),
a GPU-accelerated Optimization Engine.

The package has two components:

 - a thin wrapper around the complete C API
 - an interface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)

## Affiliation

This wrapper is developed and maintained by NVIDIA with help from the JuMP community.

## Getting Help

For assistance, please post your questions on the [JuMP community forum](https://jump.dev/forum).

If you encounter a reproducible bug, feel free to open a [GitHub issue](https://github.com/jump-dev/cuOpt.jl/issues).

## License

`cuOpt.jl` is licensed under the [Apache License, Version 2.0](https://github.com/jump-dev/cuOpt.jl/blob/master/LICENSE).

## System Requirements

Please refer to [NVIDIA cuOpt system requirements](https://docs.nvidia.com/cuopt/user-guide/latest/system-requirements.html).

## Intallation

To use cuOpt.jl, you must first separately install cuOpt.

**Installing cuOpt requires Linux.**

Please refer to the [NVIDIA cuOpt documentation](https://docs.nvidia.com/cuopt/user-guide/latest/cuopt-c/quick-start.html#installation) for installation instructions.

Please ensure the library path for `libcuopt.so` is added to `LD_LIBRARY_PATH`.

Once cuOpt is installed, add cuOpt.jl as follows:
```julia
import Pkg
Pkg.add("cuOpt")
```

## Use with JuMP

To use NVIDIA cuOpt with JuMP, use `cuOpt.Optimizer`:

```julia
using JuMP, cuOpt
model = Model(cuOpt.Optimizer)
@variable(model, x >= 0)
@variable(model, 0 <= y <= 3)
@objective(model, Min, 12x + 20y)
@constraint(model, 6x + 8y >= 100)
@constraint(model, 7x + 12y >= 120)
optimize!(model)
```

## Resources

- [NVIDIA cuOpt documentation](https://docs.nvidia.com/cuopt/user-guide/latest/)
- Julia examples can be found in the [cuOpt examples GitHub repository](https://github.com/NVIDIA/cuopt-examples/).
