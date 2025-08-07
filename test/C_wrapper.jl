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

module TestCCuOpt

using cuOpt
using Test

function small_test_model()
    cc = [1.0, -2.0]
    cl = [0.0, 0.0]
    cu = [10.0, 10.0]
    ct = [cuOpt.CUOPT_CONTINUOUS, cuOpt.CUOPT_CONTINUOUS]
    ru = [2.0, 1.0]
    rl = [0.0, 0.0]
    astart = [0, 2, 4]
    aindex = [0, 1, 0, 1]
    avalue = [1.0, 2.0, 1.0, 3.0]
    return (cc, cl, cu, ct, rl, ru, astart, aindex, avalue)
end

function test_Direct_C_call()
    (
        colcost,
        collower,
        colupper,
        coltype,
        rowlower,
        rowupper,
        astart,
        aindex,
        avalue,
    ) = small_test_model()
    n_col = convert(Cint, size(colcost, 1))
    n_row = convert(Cint, size(rowlower, 1))
    n_nz = convert(Cint, size(aindex, 1))

    colcost = convert(Array{Cdouble}, colcost)
    collower = convert(Array{Cdouble}, collower)
    colupper = convert(Array{Cdouble}, colupper)
    coltype = convert(Array{Cchar}, coltype)

    rowlower = convert(Array{Cdouble}, rowlower)
    rowupper = convert(Array{Cdouble}, rowupper)
    matstart = convert(Array{Cint}, astart)
    matindex = convert(Array{Cint}, aindex)
    matvalue = convert(Array{Cdouble}, avalue)

    colvalue, coldual =
        (Array{Cdouble,1}(undef, n_col), Array{Cdouble,1}(undef, n_col))
    rowvalue, rowdual =
        (Array{Cdouble,1}(undef, n_row), Array{Cdouble,1}(undef, n_row))
    colbasisstatus, rowbasisstatus =
        (Array{Cint,1}(undef, n_col), Array{Cint,1}(undef, n_row))

    objective_sense = Int32(cuOpt.CUOPT_MAXIMIZE)
    problem_ref = Ref{cuOpt.cuOptOptimizationProblem}()
    ret = cuOpt.cuOptCreateRangedProblem(
        n_row,
        n_col,
        objective_sense,
        0.0,
        colcost,
        matstart,
        matindex,
        matvalue,
        rowlower,
        rowupper,
        collower,
        colupper,
        coltype,
        problem_ref,
    )
    @test ret == 0
    problem = problem_ref[]

    settings_ref = Ref{cuOpt.cuOptSolverSettings}()
    ret = cuOpt.cuOptCreateSolverSettings(settings_ref)
    @test ret == 0
    settings = settings_ref[]

    solution_ref = Ref{cuOpt.cuOptSolution}()
    ret = cuOpt.cuOptSolve(problem, settings, solution_ref)
    @test ret == 0
    solution = solution_ref[]

    termination_status = Ref{Cint}(0)
    ret = cuOpt.cuOptGetTerminationStatus(solution, termination_status)
    @test ret == 0
    @test termination_status[] == cuOpt.CUOPT_TERIMINATION_STATUS_OPTIMAL

    cuOpt.cuOptDestroySolution(solution_ref)
    cuOpt.cuOptDestroySolverSettings(settings_ref)
    return cuOpt.cuOptDestroyProblem(problem_ref)
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith(string(name), "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end

TestCCuOpt.runtests()
