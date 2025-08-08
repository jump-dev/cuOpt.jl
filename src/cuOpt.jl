# SPDX-FileCopyrightText: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

module cuOpt

using Libdl

# A non-constant global. This is okay because the value gets cached in each
# ccall on the first call.
#
# We don't use `find_library` in the global scope because that path will get
# cached into the precompile file, and we won't search for libcuopt on
# subsequent `using cuOpt`. The value gets set in `__init__` instead.
global libcuopt

function __init__()
    if get(ENV, "JULIA_REGISTRYCI_AUTOMERGE", "false") == "true"
        return  # Skip. The package won't work, but it can be loaded.
    end
    libname = Libdl.find_library("libcuopt")
    if isempty(libname)
        error(
            "Could not find cuOpt library. Please ensure it is installed and in your system's library path.",
        )
    end
    global libcuopt = Libdl.dlpath(libname)
    return
end

# Handle INFINITY from C header
const INFINITY = Inf

include("gen/libcuopt.jl")
include("MOI_wrapper.jl")

import PrecompileTools

function _precompile()
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.instantiate(cuOpt.Optimizer; with_bridge_type = Float64),
    )
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 3)
    MOI.supports(model, MOI.VariableName(), typeof(x[1]))
    MOI.set(model, MOI.VariableName(), x[1], "x1")
    MOI.set(model, MOI.VariablePrimalStart(), x[1], 0.0)
    MOI.add_constraint(model, x[1], MOI.ZeroOne())
    MOI.add_constraint(model, x[2], MOI.Integer())
    for F in (MOI.VariableIndex, MOI.ScalarAffineFunction{Float64})
        MOI.supports_constraint(model, F, MOI.GreaterThan{Float64})
        MOI.supports_constraint(model, F, MOI.LessThan{Float64})
        MOI.supports_constraint(model, F, MOI.EqualTo{Float64})
    end
    MOI.supports_constraint(model, MOI.VariableIndex, MOI.ZeroOne)
    MOI.supports_constraint(model, MOI.VariableIndex, MOI.Integer)
    MOI.add_constraint(model, x[1], MOI.GreaterThan(0.0))
    MOI.add_constraint(model, x[2], MOI.LessThan(0.0))
    MOI.add_constraint(model, x[3], MOI.EqualTo(0.0))
    MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    MOI.add_constrained_variable(model, MOI.LessThan(0.0))
    MOI.add_constrained_variable(model, MOI.EqualTo(0.0))
    MOI.add_constrained_variable(model, MOI.Integer())
    MOI.add_constrained_variable(model, MOI.ZeroOne())
    set = (MOI.GreaterThan(0.0), MOI.LessThan(0.0))
    MOI.supports_add_constrained_variable(model, typeof(set))
    MOI.add_constrained_variable(model, set)
    f = 1.0 * x[1] + x[2] + x[3]
    c1 = MOI.add_constraint(model, f, MOI.GreaterThan(0.0))
    MOI.set(model, MOI.ConstraintName(), c1, "c1")
    MOI.supports(model, MOI.ConstraintName(), typeof(c1))
    MOI.add_constraint(model, f, MOI.LessThan(0.0))
    MOI.add_constraint(model, f, MOI.EqualTo(0.0))
    y, _ = MOI.add_constrained_variables(model, MOI.Nonnegatives(2))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.supports(model, MOI.ObjectiveFunction{typeof(f)}())
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    MOI.get(model, MOI.TerminationStatus())
    MOI.get(model, MOI.PrimalStatus())
    MOI.get(model, MOI.DualStatus())
    MOI.get(model, MOI.VariablePrimal(), x)
    return
end

PrecompileTools.@setup_workload begin
    PrecompileTools.@compile_workload begin
        if get(ENV, "JULIA_REGISTRYCI_AUTOMERGE", "false") != "true"
            # Try to initialize libcuopt for precompilation
            if !isdefined(cuOpt, :libcuopt)
                try
                    __init__()
                catch # If __init__ fails, skip precompilation
                    return
                end
            end
            _precompile()
        end
    end
end

end # module cuOpt
