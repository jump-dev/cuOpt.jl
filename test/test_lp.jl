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

module TestLPcuOpt

using JuMP, cuOpt
using MathOptInterface
using Test

function test_lp()
    # Create a simple linear programming model
    model = Model(cuOpt.Optimizer)

    @variable(model, x >= 0)
    @variable(model, y >= 0.5)
    @variable(model, z >= 0.5)
    @constraint(model, x + y <= 1)
    @constraint(model, x + z <= 1)
    @objective(model, Min, 2*x + y + 5)

    # Solve the model
    optimize!(model)

    # Test that the model was solved successfully
    @test termination_status(model) == MOI.OPTIMAL

    # Test objective value (expected: 5.5)
    # At optimality: x = 0, y = 0.5, z = 0.5
    # Objective: 2*0 + 0.5 + 5 = 5.5
    @test isapprox(objective_value(model), 5.5, rtol=1e-6)

    # Test variable values
    @test isapprox(value(x), 0.0, rtol=1e-6)
    @test isapprox(value(y), 0.5, rtol=1e-6)
    @test isapprox(value(z), 0.5, rtol=1e-6)

    # Test that constraints are satisfied
    @test isapprox(value(x) + value(y), 0.5, rtol=1e-6) <= 1.0
    @test isapprox(value(x) + value(z), 0.5, rtol=1e-6) <= 1.0
end

function runtests()
    test_lp()
end

end

TestLPcuOpt.runtests()
