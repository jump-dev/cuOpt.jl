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

function _test_runtests()
    model = cuOpt.Optimizer()
    MOI.Test.runtests(
        model,
        MOI.Test.Config(; atol = 1e-3, rtol = 1e-3);
        exclude = ["test_model_ScalarFunctionConstantNotZero"],
    )
    return
end

function _test_runtests_cache_optimizer()
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
            # Upstream bug: https://github.com/NVIDIA/cuopt/issues/260
            "test_constraint_ZeroOne_bounds_3",
            # Upstream bug: https://github.com/NVIDIA/cuopt/issues/112
            "test_solve_TerminationStatus_DUAL_INFEASIBLE",
            # Upstream bug: https://github.com/NVIDIA/cuopt/issues/759
            # (cuOpt crashes when given a QP with no linear constraints)
            "test_objective_qp_ObjectiveFunction_zero_ofdiag",
            "test_objective_qp_ObjectiveFunction_edge_cases",
        ],
    )
    return
end

function _test_air05()
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

function test_qp_objective()
    model = MOI.instantiate(cuOpt.Optimizer; with_cache_type = Float64)
    MOI.set(model, MOI.Silent(), true)
    x, _ = MOI.add_constrained_variable(model, (MOI.GreaterThan(0.0), MOI.LessThan(1.0)))
    y, _ = MOI.add_constrained_variable(model, (MOI.GreaterThan(0.0), MOI.LessThan(1.0)))

    # x + y == 1
    MOI.add_constraint(model,
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x), MOI.ScalarAffineTerm(1.0, y)], 0.0),
        MOI.EqualTo(1.0)
    )

    F = MOI.ScalarQuadraticFunction{Float64}
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    # 1. Homogeneous QP with diagonal objective
    # Min ¹/₂ * (x² + y²) s.t. x+y == 1
    fobj = MOI.ScalarQuadraticFunction(
        [
            MOI.ScalarQuadraticTerm(1.0, x, x),
            MOI.ScalarQuadraticTerm(1.0, y, y),
        ],
        MOI.ScalarAffineTerm{Float64}[],
        0.0,
    )
    MOI.set(model, MOI.ObjectiveFunction{F}(), fobj)
    MOI.optimize!(model)

    @test isapprox(MOI.get(model, MOI.VariablePrimal(), x), 0.5; atol=1e-4, rtol = 1e-4)
    @test isapprox(MOI.get(model, MOI.VariablePrimal(), y), 0.5; atol=1e-4, rtol = 1e-4)
    # ¹/₂ * (2 * 0.5² + 2*0.5²) == 0.5
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), 0.25; atol=1e-4, rtol = 1e-4)

    # Change diagonal coefficients
    # Min ¹/₂ * (2x² + 2y²) s.t. x+y == 1
    fobj = MOI.ScalarQuadraticFunction(
        [
            MOI.ScalarQuadraticTerm(2.0, x, x),
            MOI.ScalarQuadraticTerm(2.0, y, y),
        ],
        MOI.ScalarAffineTerm{Float64}[],
        0.0,
    )
    MOI.set(model, MOI.ObjectiveFunction{F}(), fobj)
    MOI.optimize!(model)
    # Same solution, but different objective value
    # ¹/₂ * (2 * 0.5² + 2*0.5²) == 0.5
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), 0.5; atol=1e-4, rtol = 1e-4)
    @test isapprox(MOI.get(model, MOI.VariablePrimal(), x), 0.5; atol=1e-4, rtol = 1e-4)
    @test isapprox(MOI.get(model, MOI.VariablePrimal(), y), 0.5; atol=1e-4, rtol = 1e-4)

    # QP with off-diagonal term
    # (x - y - 2)² = x² - 2xy + y² - 4x + 4y + 4
    #                   = ¹/₂ (2x² - xy - yx + 2y²) + (-4x + 4y) + 4
    # Solution is (x, y) = (1, 0)
    fobj = MOI.ScalarQuadraticFunction(
        [
            MOI.ScalarQuadraticTerm(2.0, x, x),
            MOI.ScalarQuadraticTerm(2.0, y, y),
            MOI.ScalarQuadraticTerm(-1.0, x, y),
        ],
        [
            MOI.ScalarAffineTerm(-4.0, x),
            MOI.ScalarAffineTerm(+4.0, y),
        ],
        4.0,
    )
    MOI.set(model, MOI.ObjectiveFunction{F}(), fobj)
    MOI.optimize!(model)
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), 1; atol=1e-4, rtol = 1e-4)
    @test isapprox(MOI.get(model, MOI.VariablePrimal(), x), 1.0; atol=1e-4, rtol = 1e-4)
    @test isapprox(MOI.get(model, MOI.VariablePrimal(), y), 0.0; atol=1e-4, rtol = 1e-4)

    return nothing
end

end  # TestMOIWrapper

TestMOIWrapper.runtests()
