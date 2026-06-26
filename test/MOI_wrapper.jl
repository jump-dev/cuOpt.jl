# SPDX-FileCopyrightText: Copyright (c) 2025-2026, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
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

import cuOpt
import MathOptInterface as MOI

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$name" begin
                getfield(@__MODULE__, name)()
            end
        end
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
        # cuOpt's QP solver assumes convexity, so nonconvex QCQP problems are
        # out of scope.
        exclude = ["test_quadratic_nonconvex_constraint_integration"],
    )
    return
end

function test_air05()
    src = MOI.FileFormats.MPS.Model()
    MOI.read_from_file(src, joinpath(@__DIR__, "datasets", "air05.mps"))
    model = cuOpt.Optimizer()
    MOI.set(model, MOI.RawOptimizerAttribute(cuOpt.CUOPT_LOG_TO_CONSOLE), false)
    MOI.set(model, MOI.RawOptimizerAttribute(cuOpt.CUOPT_TIME_LIMIT), 60.0)
    MOI.copy_to(model, src)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), 26374.0; rtol = 1e-4)
    return
end

# Exercise the SOC and RSOC conic conformance tests through MOI's bridge layer.
# With with_bridge_type, MOI.Bridges.full_bridge_optimizer picks up the
# SOCtoNonConvexQuadBridge / RSOCtoNonConvexQuadBridge bridges that the
# wrapper declares via ListOfNonstandardBridges. Both cones lower to
# ScalarQuadraticFunction-in-LessThan (plus an explicit t >= 0 row from the
# bridge) before reaching cuOpt, which detects the SOC structure in the
# quadratic form.
function test_runtests_bridge_optimizer()
    model = MOI.instantiate(
        cuOpt.Optimizer;
        with_cache_type = Float64,
        with_bridge_type = Float64,
    )
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
        # cuOpt's QP solver assumes convexity, so nonconvex QCQP problems are
        # out of scope.
        # Known cuOpt limitation: the barrier QCQP solver does not certify
        # primal-infeasibility or dual-infeasibility (unboundedness). It
        # returns INFEASIBLE_OR_UNBOUNDED (or NUMERICAL_ERROR in degenerate
        # cases) where MOI tests assert INFEASIBLE / DUAL_INFEASIBLE.
        # The following tests exercise this discrimination and stay excluded
        # until the cuOpt-side issue is resolved.
        #
        # The four SOC tests below all call MOI.delete on a bridged SOC
        # constraint. The SOCtoNonConvexQuadBridge cleanup triggers
        # call_in_context in MOI's variable bridge map, which crashes Julia's
        # JIT compiler (SIGSEGV in jl_type_infer) on Julia 1.12.x. Tracked
        # in https://github.com/NVIDIA/cuopt/issues/1485.
        exclude = [
            "test_quadratic_nonconvex_constraint_integration",
            r"^test_conic_RotatedSecondOrderCone_INFEASIBLE$",
            "test_conic_SecondOrderCone_no_initial_bound",
            "test_conic_SecondOrderCone_negative_initial_bound",
            "test_conic_SecondOrderCone_negative_post_bound_2",
            "test_conic_SecondOrderCone_negative_post_bound_3",
        ],
        include = [
            "test_conic_SecondOrderCone",
            "test_conic_RotatedSecondOrderCone",
        ],
    )
    return
end

end  # TestMOIWrapper

TestMOIWrapper.runtests()
