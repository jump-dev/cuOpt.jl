# How to generate the C wrapper for cuOpt

This document describes how to generate the C wrapper for cuOpt, which is for active development.

Please refer to the [main README.md](../README.md#install-nvidia-cuopt) for more information on how to install NVIDIA cuOpt.

## Building the C wrapper

Set the environment variable `CUOPT_INCLUDE_PATH` to the path of the cuOpt include directory.

You can search for the include path using the following command:

```bash
find / -path "*/cuopt/linear_programming/cuopt_c.h"
```

```bash
export CUOPT_INCLUDE_PATH=/path/to/cuopt/include
```

If the file path is something like this, ``/home/cuopt/.local/lib/python3.12/site-packages/libcuopt/include/cuopt/linear_programming/cuopt_c.h``, the path you would be setting is as followsm


```bash
export CUOPT_INCLUDE_PATH=/home/cuopt/.local/lib/python3.12/site-packages/libcuopt/include/
```

Run the following command to generate the C wrapper for cuOpt:

```bash
julia gen/gen.jl
```