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

# import Pkg; Pkg.rm("cuOpt"); Pkg.develop(path="git-repos/cuOpt.jl")
using JuMP, cuOpt

# model = Model(optimizer_with_attributes(cuOpt.Optimizer, "log_level" => 0))

model = Model(cuOpt.Optimizer)

@variable(model, x >= 0)
@variable(model, y >= 0.5)
@variable(model, z >= 0.5)
@constraint(model, x + y <= 1)
@constraint(model, x + z <= 1)
@objective(model, Min, 2*x + y + 5)

optimize!(model)

println("Objective value: ", objective_value(model))
println("x: ", value(x))
println("y: ", value(y))
println("z: ", value(z))