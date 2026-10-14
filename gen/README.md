# How to generate the C wrapper for cuOpt

This document describes how to generate the C wrapper for cuOpt, which is for active development.

Please refer to the [main README.md](../README.md#installation) for more information on how to install NVIDIA cuOpt.

## Building the C wrapper

Set the environment variable `CUOPT_INCLUDE_PATH` to the path of the cuOpt 26.10
include directory.

You can search for the include path using the following command:

```bash
find / -path "*/cuopt/mathematical_optimization/cuopt_c.h"
```

```bash
export CUOPT_INCLUDE_PATH=/path/to/cuopt/include
```

For a cuOpt source checkout, use its `cpp/include` directory, which contains
the C API header and the shared `constants.h`, `status_codes.h`, and `export.hpp`
headers:

```bash
export CUOPT_INCLUDE_PATH=/path/to/cuopt/cpp/include
```

Run the following command to generate the C wrapper for cuOpt:

```bash
julia --project=gen -e 'using Pkg; Pkg.instantiate()'
julia --project=gen gen/gen.jl
```
