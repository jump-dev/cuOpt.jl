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

module TestMOIWrapper

using Test

using cuOpt
import MathOptInterface as MOI

function runtests()
    for name in names(@__MODULE__; all = true)
        if !startswith("$(name)", "test_")
            continue
        end
        getfield(@__MODULE__, name)()
    end
    return
end

function test_runtests()
    model = cuOpt.Optimizer()

    MOI.Test.runtests(
        model,
        MOI.Test.Config(; atol = 1e-3, rtol = 1e-3);
        exclude = ["test_model_ScalarFunctionConstantNotZero"],
    )

    return
end

function test_runtests_cache_optimizer()
    model = MOI.instantiate(cuOpt.Optimizer; with_cache_type = Float64)

    MOI.Test.runtests(
        model,
        MOI.Test.Config(;
            atol = 1e-3,
            rtol = 1e-3,
            exclude = Any[
                MOI.DualObjectiveValue,
                MOI.ConstraintPrimal,
                MOI.ConstraintDual,
                MOI.ConstraintBasisStatus,
            ],
        );
        exclude = [
            "test_constraint_ZeroOne_bounds_3",
            "test_solve_TerminationStatus_DUAL_INFEASIBLE",
            "test_unbounded_MAX_SENSE",
            "test_unbounded_MIN_SENSE",
        ],
    )
    return
end

end

TestMOIWrapper.runtests()
