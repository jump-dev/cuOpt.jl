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

# Disable JuliaFormatter for this file.
#!format: off

"""
A `[`cuOptOptimizationProblem`](@ref)` object contains a representation of an LP or MIP. It is created by `[`cuOptCreateProblem`](@ref)` or `[`cuOptCreateRangedProblem`](@ref)`. It is passed to `[`cuOptSolve`](@ref)`. It should be destroyed using `[`cuOptDestroyProblem`](@ref)`.
"""
const cuOptOptimizationProblem = Ptr{Cvoid}

"""
A `[`cuOptSolverSettings`](@ref)` object contains parameter settings and other information for an LP or MIP solve. It is created by `[`cuOptCreateSolverSettings`](@ref)`. It is passed to `[`cuOptSolve`](@ref)`. It should be destroyed using `[`cuOptDestroySolverSettings`](@ref)`.
"""
const cuOptSolverSettings = Ptr{Cvoid}

"""
A `[`cuOptSolution`](@ref)` object contains the solution to an LP or MIP. It is created by `[`cuOptSolve`](@ref)`. It should be destroyed using `[`cuOptDestroySolution`](@ref)`.
"""
const cuOptSolution = Ptr{Cvoid}

"""
The type of the floating point number used by the solver. Use `[`cuOptGetFloatSize`](@ref)` to get the size of the floating point type.
"""
const cuopt_float_t = Cdouble

"""
The type of the integer number used by the solver. Use `[`cuOptGetIntSize`](@ref)` to get the size of the integer type.
"""
const cuopt_int_t = Int32

# no prototype is found for this function at cuopt_c.h:89:8, please use with caution
"""
    cuOptGetFloatSize()

Get the size of the float type.

# Returns
The size in bytes of the float type.
"""
function cuOptGetFloatSize()
    ccall((:cuOptGetFloatSize, libcuopt), Int8, ())
end

# no prototype is found for this function at cuopt_c.h:94:8, please use with caution
"""
    cuOptGetIntSize()

Get the size of the integer type used by the library.

# Returns
The size of the integer type in bytes.
"""
function cuOptGetIntSize()
    ccall((:cuOptGetIntSize, libcuopt), Int8, ())
end

"""
    cuOptGetVersion(version_major, version_minor, version_patch)

Get the version of the library.

# Arguments
* `version_major`:\\[out\\] - A pointer to a [`cuopt_int_t`](@ref) that will contain the major version number.
* `version_minor`:\\[out\\] - A pointer to a [`cuopt_int_t`](@ref) that will contain the minor version number.
* `version_patch`:\\[out\\] - A pointer to a [`cuopt_int_t`](@ref) that will contain the patch version number.
# Returns
A status code indicating success or failure.
"""
function cuOptGetVersion(version_major, version_minor, version_patch)
    ccall((:cuOptGetVersion, libcuopt), cuopt_int_t, (Ptr{cuopt_int_t}, Ptr{cuopt_int_t}, Ptr{cuopt_int_t}), version_major, version_minor, version_patch)
end

"""
    cuOptReadProblem(filename, problem_ptr)

Read an optimization problem from an MPS file.

# Arguments
* `filename`:\\[in\\] - The path to the MPS file.
* `problem_ptr`:\\[out\\] - A pointer to a [`cuOptOptimizationProblem`](@ref). On output the problem will be created and initialized with the data from the MPS file
# Returns
A status code indicating success or failure.
"""
function cuOptReadProblem(filename, problem_ptr)
    ccall((:cuOptReadProblem, libcuopt), cuopt_int_t, (Ptr{Cchar}, Ptr{cuOptOptimizationProblem}), filename, problem_ptr)
end

"""
    cuOptCreateProblem(num_constraints, num_variables, objective_sense, objective_offset, objective_coefficients, constraint_matrix_row_offsets, constraint_matrix_column_indices, constraint_matrix_coefficent_values, constraint_sense, rhs, lower_bounds, upper_bounds, variable_types, problem_ptr)

Create an optimization problem of the form

```c++
                minimize/maximize  cᵀx + offset
                  subject to       A x {=, ≤, ≥} b
                                   l ≤ x ≤ u
                                   x_i integer for some i
```

# Arguments
* `num_constraints`:\\[in\\] The number of constraints
* `num_variables`:\\[in\\] The number of variables
* `objective_sense`:\\[in\\] The objective sense ([`CUOPT_MINIMIZE`](@ref) for minimization or [`CUOPT_MAXIMIZE`](@ref) for maximization)
* `objective_offset`:\\[in\\] An offset to add to the linear objective
* `objective_coefficients`:\\[in\\] A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_variables containing the coefficients of the linear objective
* `constraint_matrix_row_offsets`:\\[in\\] A pointer to an array of type [`cuopt_int_t`](@ref) of size num\\_constraints + 1. constraint\\_matrix\\_row\\_offsets[i] is the index of the first non-zero element of the i-th constraint in constraint\\_matrix\\_column\\_indices and constraint\\_matrix\\_coefficent\\_values. This is part of the compressed sparse row representation of the constraint matrix
* `constraint_matrix_column_indices`:\\[in\\] A pointer to an array of type [`cuopt_int_t`](@ref) of size constraint\\_matrix\\_row\\_offsets[num\\_constraints] containing the column indices of the non-zero elements of the constraint matrix. This is part of the compressed sparse row representation of the constraint matrix
* `constraint_matrix_coefficent_values`:\\[in\\] A pointer to an array of type [`cuopt_float_t`](@ref) of size constraint\\_matrix\\_row\\_offsets[num\\_constraints] containing the values of the non-zero elements of the constraint matrix. This is part of the compressed sparse row representation of the constraint matrix
* `constraint_sense`:\\[in\\] A pointer to an array of type char of size num\\_constraints containing the sense of the constraints ([`CUOPT_LESS_THAN`](@ref), [`CUOPT_GREATER_THAN`](@ref), or [`CUOPT_EQUAL`](@ref))
* `rhs`:\\[in\\] A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_constraints containing the right-hand side of the constraints
* `lower_bounds`:\\[in\\] A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_variables containing the lower bounds of the variables
* `upper_bounds`:\\[in\\] A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_variables containing the upper bounds of the variables
* `variable_types`:\\[in\\] A pointer to an array of type char of size num\\_variables containing the types of the variables ([`CUOPT_CONTINUOUS`](@ref) or [`CUOPT_INTEGER`](@ref))
* `problem_ptr`:\\[out\\] Pointer to store the created optimization problem
# Returns
[`CUOPT_SUCCESS`](@ref) if successful, CUOPT\\_ERROR otherwise
"""
function cuOptCreateProblem(num_constraints, num_variables, objective_sense, objective_offset, objective_coefficients, constraint_matrix_row_offsets, constraint_matrix_column_indices, constraint_matrix_coefficent_values, constraint_sense, rhs, lower_bounds, upper_bounds, variable_types, problem_ptr)
    ccall((:cuOptCreateProblem, libcuopt), cuopt_int_t, (cuopt_int_t, cuopt_int_t, cuopt_int_t, cuopt_float_t, Ptr{cuopt_float_t}, Ptr{cuopt_int_t}, Ptr{cuopt_int_t}, Ptr{cuopt_float_t}, Ptr{Cchar}, Ptr{cuopt_float_t}, Ptr{cuopt_float_t}, Ptr{cuopt_float_t}, Ptr{Cchar}, Ptr{cuOptOptimizationProblem}), num_constraints, num_variables, objective_sense, objective_offset, objective_coefficients, constraint_matrix_row_offsets, constraint_matrix_column_indices, constraint_matrix_coefficent_values, constraint_sense, rhs, lower_bounds, upper_bounds, variable_types, problem_ptr)
end

"""
    cuOptCreateRangedProblem(num_constraints, num_variables, objective_sense, objective_offset, objective_coefficients, constraint_matrix_row_offsets, constraint_matrix_column_indices, constraint_matrix_coefficients, constraint_lower_bounds, constraint_upper_bounds, variable_lower_bounds, variable_upper_bounds, variable_types, problem_ptr)

Create an optimization problem of the form *

```c++
                minimize/maximize  c^T x + offset
                  subject to       bl <= A*x <= bu
                                   l <= x <= u
                                   x_i integer for some i
```

# Arguments
* `num_constraints`:\\[in\\] - The number of constraints.
* `num_variables`:\\[in\\] - The number of variables.
* `objective_sense`:\\[in\\] - The objective sense ([`CUOPT_MINIMIZE`](@ref) for minimization or [`CUOPT_MAXIMIZE`](@ref) for maximization)
* `objective_offset`:\\[in\\] - An offset to add to the linear objective.
* `objective_coefficients`:\\[in\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_variables containing the coefficients of the linear objective.
* `constraint_matrix_row_offsets`:\\[in\\] - A pointer to an array of type [`cuopt_int_t`](@ref) of size num\\_constraints + 1. constraint\\_matrix\\_row\\_offsets[i] is the index of the first non-zero element of the i-th constraint in constraint\\_matrix\\_column\\_indices and constraint\\_matrix\\_coefficients.
* `constraint_matrix_column_indices`:\\[in\\] - A pointer to an array of type [`cuopt_int_t`](@ref) of size constraint\\_matrix\\_row\\_offsets[num\\_constraints] containing the column indices of the non-zero elements of the constraint matrix.
* `constraint_matrix_coefficients`:\\[in\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size constraint\\_matrix\\_row\\_offsets[num\\_constraints] containing the values of the non-zero elements of the constraint matrix.
* `constraint_lower_bounds`:\\[in\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_constraints containing the lower bounds of the constraints.
* `constraint_upper_bounds`:\\[in\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_constraints containing the upper bounds of the constraints.
* `variable_lower_bounds`:\\[in\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_variables containing the lower bounds of the variables.
* `variable_upper_bounds`:\\[in\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_variables containing the upper bounds of the variables.
* `variable_types`:\\[in\\] - A pointer to an array of type char of size num\\_variables containing the types of the variables ([`CUOPT_CONTINUOUS`](@ref) or [`CUOPT_INTEGER`](@ref)).
* `problem_ptr`:\\[out\\] - A pointer to a [`cuOptOptimizationProblem`](@ref). On output the problem will be created and initialized with the provided data.
# Returns
A status code indicating success or failure.
"""
function cuOptCreateRangedProblem(num_constraints, num_variables, objective_sense, objective_offset, objective_coefficients, constraint_matrix_row_offsets, constraint_matrix_column_indices, constraint_matrix_coefficients, constraint_lower_bounds, constraint_upper_bounds, variable_lower_bounds, variable_upper_bounds, variable_types, problem_ptr)
    ccall((:cuOptCreateRangedProblem, libcuopt), cuopt_int_t, (cuopt_int_t, cuopt_int_t, cuopt_int_t, cuopt_float_t, Ptr{cuopt_float_t}, Ptr{cuopt_int_t}, Ptr{cuopt_int_t}, Ptr{cuopt_float_t}, Ptr{cuopt_float_t}, Ptr{cuopt_float_t}, Ptr{cuopt_float_t}, Ptr{cuopt_float_t}, Ptr{Cchar}, Ptr{cuOptOptimizationProblem}), num_constraints, num_variables, objective_sense, objective_offset, objective_coefficients, constraint_matrix_row_offsets, constraint_matrix_column_indices, constraint_matrix_coefficients, constraint_lower_bounds, constraint_upper_bounds, variable_lower_bounds, variable_upper_bounds, variable_types, problem_ptr)
end

"""
    cuOptDestroyProblem(problem_ptr)

Destroy an optimization problem

# Arguments
* `problem_ptr`:\\[in,out\\] - A pointer to a [`cuOptOptimizationProblem`](@ref). On output the problem will be destroyed, and the pointer will be set to NULL.
"""
function cuOptDestroyProblem(problem_ptr)
    ccall((:cuOptDestroyProblem, libcuopt), Cvoid, (Ptr{cuOptOptimizationProblem},), problem_ptr)
end

"""
    cuOptGetNumConstraints(problem, num_constraints_ptr)

Get the number of constraints of an optimization problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `num_constraints_ptr`:\\[out\\] - A pointer to a [`cuopt_int_t`](@ref) that will contain the number of constraints on output.
# Returns
A status code indicating success or failure.
"""
function cuOptGetNumConstraints(problem, num_constraints_ptr)
    ccall((:cuOptGetNumConstraints, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{cuopt_int_t}), problem, num_constraints_ptr)
end

"""
    cuOptGetNumVariables(problem, num_variables_ptr)

Get the number of variables of an optimization problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `num_variables_ptr`:\\[out\\] - A pointer to a [`cuopt_int_t`](@ref) that will contain the number of variables on output.
# Returns
A status code indicating success or failure.
"""
function cuOptGetNumVariables(problem, num_variables_ptr)
    ccall((:cuOptGetNumVariables, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{cuopt_int_t}), problem, num_variables_ptr)
end

"""
    cuOptGetObjectiveSense(problem, objective_sense_ptr)

Get the objective sense of an optimization problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `objective_sense_ptr`:\\[out\\] - A pointer to a [`cuopt_int_t`](@ref) that on output will contain the objective sense.
# Returns
A status code indicating success or failure.
"""
function cuOptGetObjectiveSense(problem, objective_sense_ptr)
    ccall((:cuOptGetObjectiveSense, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{cuopt_int_t}), problem, objective_sense_ptr)
end

"""
    cuOptGetObjectiveOffset(problem, objective_offset_ptr)

Get the objective offset of an optimization problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `objective_offset_ptr`:\\[out\\] - A pointer to a [`cuopt_float_t`](@ref) that on output will contain the objective offset.
# Returns
A status code indicating success or failure.
"""
function cuOptGetObjectiveOffset(problem, objective_offset_ptr)
    ccall((:cuOptGetObjectiveOffset, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{cuopt_float_t}), problem, objective_offset_ptr)
end

"""
    cuOptGetObjectiveCoefficients(problem, objective_coefficients_ptr)

Get the objective coefficients of an optimization problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `objective_coefficients_ptr`:\\[out\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_variables that on output will contain the objective coefficients.
# Returns
A status code indicating success or failure.
"""
function cuOptGetObjectiveCoefficients(problem, objective_coefficients_ptr)
    ccall((:cuOptGetObjectiveCoefficients, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{cuopt_float_t}), problem, objective_coefficients_ptr)
end

"""
    cuOptGetNumNonZeros(problem, num_non_zeros_ptr)

Get the number of non-zero elements in the constraint matrix of an optimization problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `num_non_zeros_ptr`:\\[out\\] - A pointer to a [`cuopt_int_t`](@ref) that on output will contain the number of non-zeros in the constraint matrix.
# Returns
A status code indicating success or failure.
"""
function cuOptGetNumNonZeros(problem, num_non_zeros_ptr)
    ccall((:cuOptGetNumNonZeros, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{cuopt_int_t}), problem, num_non_zeros_ptr)
end

"""
    cuOptGetConstraintMatrix(problem, constraint_matrix_row_offsets_ptr, constraint_matrix_column_indices_ptr, constraint_matrix_coefficients_ptr)

Get the constraint matrix of an optimization problem in compressed sparse row format.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `constraint_matrix_row_offsets_ptr`:\\[out\\] - A pointer to an array of type [`cuopt_int_t`](@ref) of size num\\_constraints + 1 that on output will contain the row offsets of the constraint matrix.
* `constraint_matrix_column_indices_ptr`:\\[out\\] - A pointer to an array of type [`cuopt_int_t`](@ref) of size equal to the number of nonzeros that on output will contain the column indices of the non-zero entries of the constraint matrix.
* `constraint_matrix_coefficients_ptr`:\\[out\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size equal to the number of nonzeros that on output will contain the coefficients of the non-zero entries of the constraint matrix.
# Returns
A status code indicating success or failure.
"""
function cuOptGetConstraintMatrix(problem, constraint_matrix_row_offsets_ptr, constraint_matrix_column_indices_ptr, constraint_matrix_coefficients_ptr)
    ccall((:cuOptGetConstraintMatrix, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{cuopt_int_t}, Ptr{cuopt_int_t}, Ptr{cuopt_float_t}), problem, constraint_matrix_row_offsets_ptr, constraint_matrix_column_indices_ptr, constraint_matrix_coefficients_ptr)
end

"""
    cuOptGetConstraintSense(problem, constraint_sense_ptr)

Get the constraint sense of an optimization problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `constraint_sense_ptr`:\\[out\\] - A pointer to an array of type char of size num\\_constraints that on output will contain the sense of the constraints.
# Returns
A status code indicating success or failure.
"""
function cuOptGetConstraintSense(problem, constraint_sense_ptr)
    ccall((:cuOptGetConstraintSense, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{Cchar}), problem, constraint_sense_ptr)
end

"""
    cuOptGetConstraintRightHandSide(problem, rhs_ptr)

Get the right-hand side of an optimization problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `rhs_ptr`:\\[out\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_constraints that on output will contain the right-hand side of the constraints.
# Returns
A status code indicating success or failure.
"""
function cuOptGetConstraintRightHandSide(problem, rhs_ptr)
    ccall((:cuOptGetConstraintRightHandSide, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{cuopt_float_t}), problem, rhs_ptr)
end

"""
    cuOptGetConstraintLowerBounds(problem, lower_bounds_ptr)

Get the lower bounds of an optimization problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `lower_bounds_ptr`:\\[out\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_constraints that on output will contain the lower bounds of the constraints.
# Returns
A status code indicating success or failure.
"""
function cuOptGetConstraintLowerBounds(problem, lower_bounds_ptr)
    ccall((:cuOptGetConstraintLowerBounds, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{cuopt_float_t}), problem, lower_bounds_ptr)
end

"""
    cuOptGetConstraintUpperBounds(problem, upper_bounds_ptr)

Get the upper bounds of an optimization problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `upper_bounds_ptr`:\\[out\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_constraints that on output will contain the upper bounds of the constraints.
# Returns
A status code indicating success or failure.
"""
function cuOptGetConstraintUpperBounds(problem, upper_bounds_ptr)
    ccall((:cuOptGetConstraintUpperBounds, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{cuopt_float_t}), problem, upper_bounds_ptr)
end

"""
    cuOptGetVariableLowerBounds(problem, lower_bounds_ptr)

Get the lower bounds of an optimization problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `lower_bounds_ptr`:\\[out\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_variables that on output will contain the lower bounds of the variables.
# Returns
A status code indicating success or failure.
"""
function cuOptGetVariableLowerBounds(problem, lower_bounds_ptr)
    ccall((:cuOptGetVariableLowerBounds, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{cuopt_float_t}), problem, lower_bounds_ptr)
end

"""
    cuOptGetVariableUpperBounds(problem, upper_bounds_ptr)

Get the upper bounds of an optimization problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `upper_bounds_ptr`:\\[out\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_variables that on output will contain the upper bounds of the variables.
# Returns
A status code indicating success or failure.
"""
function cuOptGetVariableUpperBounds(problem, upper_bounds_ptr)
    ccall((:cuOptGetVariableUpperBounds, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{cuopt_float_t}), problem, upper_bounds_ptr)
end

"""
    cuOptGetVariableTypes(problem, variable_types_ptr)

Get the variable types of an optimization problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `variable_types_ptr`:\\[out\\] - A pointer to an array of type char of size num\\_variables that on output will contain the types of the variables ([`CUOPT_CONTINUOUS`](@ref) or [`CUOPT_INTEGER`](@ref)).
# Returns
A status code indicating success or failure.
"""
function cuOptGetVariableTypes(problem, variable_types_ptr)
    ccall((:cuOptGetVariableTypes, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{Cchar}), problem, variable_types_ptr)
end

"""
    cuOptCreateSolverSettings(settings_ptr)

Create a solver settings object.

# Arguments
* `settings_ptr`:\\[out\\] - A pointer to a [`cuOptSolverSettings`](@ref) object. On output the solver settings will be created and initialized.
# Returns
A status code indicating success or failure.
"""
function cuOptCreateSolverSettings(settings_ptr)
    ccall((:cuOptCreateSolverSettings, libcuopt), cuopt_int_t, (Ptr{cuOptSolverSettings},), settings_ptr)
end

"""
    cuOptDestroySolverSettings(settings_ptr)

Destroy a solver settings object.

# Arguments
* `settings_ptr`:\\[in,out\\] - A pointer to a [`cuOptSolverSettings`](@ref) object. On output the solver settings will be destroyed and the pointer will be set to NULL.
"""
function cuOptDestroySolverSettings(settings_ptr)
    ccall((:cuOptDestroySolverSettings, libcuopt), Cvoid, (Ptr{cuOptSolverSettings},), settings_ptr)
end

"""
    cuOptSetParameter(settings, parameter_name, parameter_value)

Set a parameter of a solver settings object.

# Arguments
* `settings`:\\[in\\] - The solver settings object.
* `parameter_name`:\\[in\\] - The name of the parameter to set.
* `parameter_value`:\\[in\\] - The value of the parameter to set.
"""
function cuOptSetParameter(settings, parameter_name, parameter_value)
    ccall((:cuOptSetParameter, libcuopt), cuopt_int_t, (cuOptSolverSettings, Ptr{Cchar}, Ptr{Cchar}), settings, parameter_name, parameter_value)
end

"""
    cuOptGetParameter(settings, parameter_name, parameter_value_size, parameter_value)

Get a parameter of a solver settings object.

# Arguments
* `settings`:\\[in\\] - The solver settings object.
* `parameter_name`:\\[in\\] - The name of the parameter to get.
* `parameter_value_size`:\\[in\\] - The size of the parameter value buffer.
* `parameter_value`:\\[out\\] - A pointer to an array of characters that on output will contain the value of the parameter.
# Returns
A status code indicating success or failure.
"""
function cuOptGetParameter(settings, parameter_name, parameter_value_size, parameter_value)
    ccall((:cuOptGetParameter, libcuopt), cuopt_int_t, (cuOptSolverSettings, Ptr{Cchar}, cuopt_int_t, Ptr{Cchar}), settings, parameter_name, parameter_value_size, parameter_value)
end

"""
    cuOptSetIntegerParameter(settings, parameter_name, parameter_value)

Set an integer parameter of a solver settings object.

# Arguments
* `settings`:\\[in\\] - The solver settings object.
* `parameter_name`:\\[in\\] - The name of the parameter to set.
* `parameter_value`:\\[in\\] - The value of the parameter to set.
# Returns
A status code indicating success or failure.
"""
function cuOptSetIntegerParameter(settings, parameter_name, parameter_value)
    ccall((:cuOptSetIntegerParameter, libcuopt), cuopt_int_t, (cuOptSolverSettings, Ptr{Cchar}, cuopt_int_t), settings, parameter_name, parameter_value)
end

"""
    cuOptGetIntegerParameter(settings, parameter_name, parameter_value)

Get an integer parameter of a solver settings object.

# Arguments
* `settings`:\\[in\\] - The solver settings object.
* `parameter_name`:\\[in\\] - The name of the parameter to get.
* `parameter_value`:\\[out\\] - A pointer to a [`cuopt_int_t`](@ref) that on output will contain the value of the parameter.
# Returns
A status code indicating success or failure.
"""
function cuOptGetIntegerParameter(settings, parameter_name, parameter_value)
    ccall((:cuOptGetIntegerParameter, libcuopt), cuopt_int_t, (cuOptSolverSettings, Ptr{Cchar}, Ptr{cuopt_int_t}), settings, parameter_name, parameter_value)
end

"""
    cuOptSetFloatParameter(settings, parameter_name, parameter_value)

Set a float parameter of a solver settings object.

# Arguments
* `settings`:\\[in\\] - The solver settings object.
* `parameter_name`:\\[in\\] - The name of the parameter to set.
* `parameter_value`:\\[in\\] - The value of the parameter to set.
# Returns
A status code indicating success or failure.
"""
function cuOptSetFloatParameter(settings, parameter_name, parameter_value)
    ccall((:cuOptSetFloatParameter, libcuopt), cuopt_int_t, (cuOptSolverSettings, Ptr{Cchar}, cuopt_float_t), settings, parameter_name, parameter_value)
end

"""
    cuOptGetFloatParameter(settings, parameter_name, parameter_value)

Get a float parameter of a solver settings object.

# Arguments
* `settings`:\\[in\\] - The solver settings object.
* `parameter_name`:\\[in\\] - The name of the parameter to get.
* `parameter_value`:\\[out\\] - A pointer to a [`cuopt_float_t`](@ref) that on output will contain the value of the parameter.
# Returns
A status code indicating success or failure.
"""
function cuOptGetFloatParameter(settings, parameter_name, parameter_value)
    ccall((:cuOptGetFloatParameter, libcuopt), cuopt_int_t, (cuOptSolverSettings, Ptr{Cchar}, Ptr{cuopt_float_t}), settings, parameter_name, parameter_value)
end

"""
    cuOptIsMIP(problem, is_mip_ptr)

Check if an optimization problem is a mixed integer programming problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `is_mip_ptr`:\\[out\\] - A pointer to a [`cuopt_int_t`](@ref) that on output will be 0 if the problem contains only continuous variables, or 1 if the problem contains integer variables.
# Returns
A status code indicating success or failure.
"""
function cuOptIsMIP(problem, is_mip_ptr)
    ccall((:cuOptIsMIP, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, Ptr{cuopt_int_t}), problem, is_mip_ptr)
end

"""
    cuOptSolve(problem, settings, solution_ptr)

Solve an optimization problem.

# Arguments
* `problem`:\\[in\\] - The optimization problem.
* `settings`:\\[in\\] - The solver settings.
* `solution_ptr`:\\[out\\] - A pointer to a [`cuOptSolution`](@ref) object. On output the solution will be created.
# Returns
A status code indicating success or failure.
"""
function cuOptSolve(problem, settings, solution_ptr)
    ccall((:cuOptSolve, libcuopt), cuopt_int_t, (cuOptOptimizationProblem, cuOptSolverSettings, Ptr{cuOptSolution}), problem, settings, solution_ptr)
end

"""
    cuOptDestroySolution(solution_ptr)

Destroy a solution object.

# Arguments
* `solution_ptr`:\\[in,out\\] - A pointer to a [`cuOptSolution`](@ref) object. On output the solution will be destroyed and the pointer will be set to NULL.
"""
function cuOptDestroySolution(solution_ptr)
    ccall((:cuOptDestroySolution, libcuopt), Cvoid, (Ptr{cuOptSolution},), solution_ptr)
end

"""
    cuOptGetTerminationStatus(solution, termination_status_ptr)

Get the termination reason of an optimization problem.

# Arguments
* `solution`:\\[in\\] - The solution object.
* `termination_reason_ptr`:\\[out\\] - A pointer to a [`cuopt_int_t`](@ref) that on output will contain the termination reason.
# Returns
A status code indicating success or failure.
"""
function cuOptGetTerminationStatus(solution, termination_status_ptr)
    ccall((:cuOptGetTerminationStatus, libcuopt), cuopt_int_t, (cuOptSolution, Ptr{cuopt_int_t}), solution, termination_status_ptr)
end

function cuOptGetErrorStatus(solution, error_status_ptr)
    ccall((:cuOptGetErrorStatus, libcuopt), cuopt_int_t, (cuOptSolution, Ptr{cuopt_int_t}), solution, error_status_ptr)
end

function cuOptGetErrorString(solution, error_string_ptr, error_string_size)
    ccall((:cuOptGetErrorString, libcuopt), cuopt_int_t, (cuOptSolution, Ptr{Cchar}, cuopt_int_t), solution, error_string_ptr, error_string_size)
end

function cuOptGetPrimalSolution(solution, solution_values)
    ccall((:cuOptGetPrimalSolution, libcuopt), cuopt_int_t, (cuOptSolution, Ptr{cuopt_float_t}), solution, solution_values)
end

"""
    cuOptGetObjectiveValue(solution, objective_value_ptr)

Get the objective value of an optimization problem.

# Arguments
* `solution`:\\[in\\] - The solution object.
* `objective_value_ptr`:\\[in,out\\] - A pointer to a [`cuopt_float_t`](@ref) that will contain the objective value.
# Returns
A status code indicating success or failure.
"""
function cuOptGetObjectiveValue(solution, objective_value_ptr)
    ccall((:cuOptGetObjectiveValue, libcuopt), cuopt_int_t, (cuOptSolution, Ptr{cuopt_float_t}), solution, objective_value_ptr)
end

"""
    cuOptGetSolveTime(solution, solve_time_ptr)

Get the solve time of an optimization problem.

# Arguments
* `solution`:\\[in\\] - The solution object.
* `solve_time_ptr`:\\[in,out\\] - A pointer to a [`cuopt_float_t`](@ref) that will contain the solve time.
# Returns
A status code indicating success or failure.
"""
function cuOptGetSolveTime(solution, solve_time_ptr)
    ccall((:cuOptGetSolveTime, libcuopt), cuopt_int_t, (cuOptSolution, Ptr{cuopt_float_t}), solution, solve_time_ptr)
end

"""
    cuOptGetMIPGap(solution, mip_gap_ptr)

Get the relative MIP gap of an optimization problem.

# Arguments
* `solution`:\\[in\\] - The solution object.
* `mip_gap_ptr`:\\[in,out\\] - A pointer to a [`cuopt_float_t`](@ref) that will contain the relative MIP gap.
# Returns
A status code indicating success or failure.
"""
function cuOptGetMIPGap(solution, mip_gap_ptr)
    ccall((:cuOptGetMIPGap, libcuopt), cuopt_int_t, (cuOptSolution, Ptr{cuopt_float_t}), solution, mip_gap_ptr)
end

"""
    cuOptGetSolutionBound(solution, solution_bound_ptr)

Get the solution bound of an optimization problem.

# Arguments
* `solution`:\\[in\\] - The solution object.
* `solution_bound_ptr`:\\[in,out\\] - A pointer to a [`cuopt_float_t`](@ref) that will contain the solution bound.
# Returns
A status code indicating success or failure.
"""
function cuOptGetSolutionBound(solution, solution_bound_ptr)
    ccall((:cuOptGetSolutionBound, libcuopt), cuopt_int_t, (cuOptSolution, Ptr{cuopt_float_t}), solution, solution_bound_ptr)
end

"""
    cuOptGetDualSolution(solution, dual_solution_ptr)

Get the dual solution of an optimization problem.

# Arguments
* `solution`:\\[in\\] - The solution object.
* `dual_solution_ptr`:\\[in,out\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_constraints that will contain the dual solution.
# Returns
A status code indicating success or failure.
"""
function cuOptGetDualSolution(solution, dual_solution_ptr)
    ccall((:cuOptGetDualSolution, libcuopt), cuopt_int_t, (cuOptSolution, Ptr{cuopt_float_t}), solution, dual_solution_ptr)
end

"""
    cuOptGetDualObjectiveValue(solution, dual_objective_value_ptr)

Get the dual objective value of an optimization problem.

# Arguments
* `solution`:\\[in\\] - The solution object.
* `dual_objective_value_ptr`:\\[in,out\\] - A pointer to a [`cuopt_float_t`](@ref) that will contain the dual objective value.
# Returns
A status code indicating success or failure.
"""
function cuOptGetDualObjectiveValue(solution, dual_objective_value_ptr)
    ccall((:cuOptGetDualObjectiveValue, libcuopt), cuopt_int_t, (cuOptSolution, Ptr{cuopt_float_t}), solution, dual_objective_value_ptr)
end

"""
    cuOptGetReducedCosts(solution, reduced_cost_ptr)

Get the reduced costs of an optimization problem.

# Arguments
* `solution`:\\[in\\] - The solution object.
* `reduced_cost_ptr`:\\[in,out\\] - A pointer to an array of type [`cuopt_float_t`](@ref) of size num\\_variables that will contain the reduced cost.
# Returns
A status code indicating success or failure.
"""
function cuOptGetReducedCosts(solution, reduced_cost_ptr)
    ccall((:cuOptGetReducedCosts, libcuopt), cuopt_int_t, (cuOptSolution, Ptr{cuopt_float_t}), solution, reduced_cost_ptr)
end

const CUOPT_INSTANTIATE_FLOAT = 0

const CUOPT_INSTANTIATE_DOUBLE = 1

const CUOPT_INSTANTIATE_INT32 = 1

const CUOPT_INSTANTIATE_INT64 = 0

const CUOPT_ABSOLUTE_DUAL_TOLERANCE = "absolute_dual_tolerance"

const CUOPT_RELATIVE_DUAL_TOLERANCE = "relative_dual_tolerance"

const CUOPT_ABSOLUTE_PRIMAL_TOLERANCE = "absolute_primal_tolerance"

const CUOPT_RELATIVE_PRIMAL_TOLERANCE = "relative_primal_tolerance"

const CUOPT_ABSOLUTE_GAP_TOLERANCE = "absolute_gap_tolerance"

const CUOPT_RELATIVE_GAP_TOLERANCE = "relative_gap_tolerance"

const CUOPT_INFEASIBILITY_DETECTION = "infeasibility_detection"

const CUOPT_STRICT_INFEASIBILITY = "strict_infeasibility"

const CUOPT_PRIMAL_INFEASIBLE_TOLERANCE = "primal_infeasible_tolerance"

const CUOPT_DUAL_INFEASIBLE_TOLERANCE = "dual_infeasible_tolerance"

const CUOPT_ITERATION_LIMIT = "iteration_limit"

const CUOPT_TIME_LIMIT = "time_limit"

const CUOPT_PDLP_SOLVER_MODE = "pdlp_solver_mode"

const CUOPT_METHOD = "method"

const CUOPT_PER_CONSTRAINT_RESIDUAL = "per_constraint_residual"

const CUOPT_SAVE_BEST_PRIMAL_SO_FAR = "save_best_primal_so_far"

const CUOPT_FIRST_PRIMAL_FEASIBLE = "first_primal_feasible"

const CUOPT_LOG_FILE = "log_file"

const CUOPT_LOG_TO_CONSOLE = "log_to_console"

const CUOPT_CROSSOVER = "crossover"

const CUOPT_FOLDING = "folding"

const CUOPT_AUGMENTED = "augmented"

const CUOPT_DUALIZE = "dualize"

const CUOPT_ORDERING = "ordering"

const CUOPT_BARRIER_DUAL_INITIAL_POINT = "barrier_dual_initial_point"

const CUOPT_ELIMINATE_DENSE_COLUMNS = "eliminate_dense_columns"

const CUOPT_CUDSS_DETERMINISTIC = "cudss_deterministic"

const CUOPT_PRESOLVE = "presolve"

const CUOPT_DUAL_POSTSOLVE = "dual_postsolve"

const CUOPT_MIP_ABSOLUTE_TOLERANCE = "mip_absolute_tolerance"

const CUOPT_MIP_RELATIVE_TOLERANCE = "mip_relative_tolerance"

const CUOPT_MIP_INTEGRALITY_TOLERANCE = "mip_integrality_tolerance"

const CUOPT_MIP_ABSOLUTE_GAP = "mip_absolute_gap"

const CUOPT_MIP_RELATIVE_GAP = "mip_relative_gap"

const CUOPT_MIP_HEURISTICS_ONLY = "mip_heuristics_only"

const CUOPT_MIP_SCALING = "mip_scaling"

const CUOPT_MIP_PRESOLVE = "mip_presolve"

const CUOPT_SOLUTION_FILE = "solution_file"

const CUOPT_NUM_CPU_THREADS = "num_cpu_threads"

const CUOPT_USER_PROBLEM_FILE = "user_problem_file"

const CUOPT_TERIMINATION_STATUS_NO_TERMINATION = 0

const CUOPT_TERIMINATION_STATUS_OPTIMAL = 1

const CUOPT_TERIMINATION_STATUS_INFEASIBLE = 2

const CUOPT_TERIMINATION_STATUS_UNBOUNDED = 3

const CUOPT_TERIMINATION_STATUS_ITERATION_LIMIT = 4

const CUOPT_TERIMINATION_STATUS_TIME_LIMIT = 5

const CUOPT_TERIMINATION_STATUS_NUMERICAL_ERROR = 6

const CUOPT_TERIMINATION_STATUS_PRIMAL_FEASIBLE = 7

const CUOPT_TERIMINATION_STATUS_FEASIBLE_FOUND = 8

const CUOPT_TERIMINATION_STATUS_CONCURRENT_LIMIT = 9

const CUOPT_MINIMIZE = 1

const CUOPT_MAXIMIZE = -1

const CUOPT_LESS_THAN = Cchar('L')

const CUOPT_GREATER_THAN = Cchar('G')

const CUOPT_EQUAL = Cchar('E')

const CUOPT_CONTINUOUS = Cchar('C')

const CUOPT_INTEGER = Cchar('I')

const CUOPT_INFINITY = INFINITY

const CUOPT_PDLP_SOLVER_MODE_STABLE1 = 0

const CUOPT_PDLP_SOLVER_MODE_STABLE2 = 1

const CUOPT_PDLP_SOLVER_MODE_METHODICAL1 = 2

const CUOPT_PDLP_SOLVER_MODE_FAST1 = 3

const CUOPT_PDLP_SOLVER_MODE_STABLE3 = 4

const CUOPT_METHOD_CONCURRENT = 0

const CUOPT_METHOD_PDLP = 1

const CUOPT_METHOD_DUAL_SIMPLEX = 2

const CUOPT_METHOD_BARRIER = 3

const CUOPT_SUCCESS = 0

const CUOPT_INVALID_ARGUMENT = 1

const CUOPT_MPS_FILE_ERROR = 2

const CUOPT_MPS_PARSE_ERROR = 3

const CUOPT_VALIDATION_ERROR = 4

const CUOPT_OUT_OF_MEMORY = 5

const CUOPT_RUNTIME_ERROR = 6
