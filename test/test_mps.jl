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

module TestMPScuOpt

using JuMP, cuOpt
using MathOptInterface
using Test
using Base.Filesystem

function test_mps()
    mps_file = "air05.mps"
    full_path = joinpath(@__DIR__, "datasets", mps_file)
    model = read_from_file(full_path)
    @test model isa JuMP.Model
    set_optimizer(model, cuOpt.Optimizer)

    set_attribute(model, cuOpt.CUOPT_LOG_FILE, "$(splitext(mps_file)[1]).log")
    set_attribute(model, cuOpt.CUOPT_LOG_TO_CONSOLE, false)
    set_attribute(model, cuOpt.CUOPT_TIME_LIMIT, 60.0)

    # solve
    optimize!(model)

    # We should get optimal solution for this problem
    @test termination_status(model) == MOI.OPTIMAL
    @test isapprox(objective_value(model), 26374.0, rtol = 1e-4)
end

function runtests()
    return test_mps()
end

end

TestMPScuOpt.runtests()
