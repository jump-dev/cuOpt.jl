# cuOpt.jl — Julia wrapper for [NVIDIA cuOpt](https://github.com/NVIDIA/cuOpt), a GPU-accelerated Optimization Engine

## License

cuOpt.jl is licensed under the Apache License, Version 2.0. See the [LICENSE](LICENSE) file for details.

## Getting Help

For assistance, please post your questions on the [JuMP community forum](https://jump.dev/forum).

If you encounter a reproducible bug, feel free to open a [GitHub issue](https://github.com/jump-dev/cuOpt.jl/issues).

## System Requirements

Please refer to [NVIDIA cuOpt system requirements](https://docs.nvidia.com/cuopt/user-guide/latest/system-requirements.html).


## Intallation

### Install NVIDIA cuOpt

Please refer to the [NVIDIA cuOpt installation documentation](https://docs.nvidia.com/cuopt/user-guide/latest/cuopt-c/quick-start.html#installation) for instructions on installing cuOpt.

This ensures all the dependencies that are erquired will be installed along with cuOpt.

Please ensure the library path for ``libcuopt.so`` is added to ``LD_LIBRARY_PATH``.


### Adding cuopt.jl package

```julia
import Pkg
Pkg.add("cuOpt")
```

### Example

To use NVIDIA cuOpt with JUMP, use ``cuOpt.Optimizer``

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

### Local installation

If you want to install cuOpt.jl from local github repo, please complete all other parts of installing NVIDIA cuOpt and adding to LD_LIBRARY_PATH.

Clone this repository to your local system and follow the steps below to use cuOpt.jl.

```julia
import Pkg
Pkg.add("JuMP")

# Add local pacakge
Pkg.develop(path="/path/to/cuOpt.jl")

using JuMP
using cuOpt

model = Model(cuOpt.Optimizer)
```

## Resources

- [NVIDIA cuOpt documentation](https://docs.nvidia.com/cuopt/user-guide/latest/)
- Julia examples can be found in the [cuOpt examples GitHub repository](https://github.com/NVIDIA/cuopt-examples/).











