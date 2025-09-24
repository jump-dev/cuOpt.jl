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

# This file is a modification of the MOI wrapper for HiGHS.jl, which can be
# found at https://github.com/jump-dev/HiGHS.jl/blob/master/src/MOI_wrapper.jl
# The HiGHS wrapper is released under an MIT license, a copy of which can be
# found in `/thirdparty/THIRD_PARTY_LICENSES` or at https://opensource.org/licenses/MIT.

import MathOptInterface as MOI
const CleverDicts = MOI.Utilities.CleverDicts

const _SCALAR_SETS =
    Union{MOI.GreaterThan{Float64},MOI.LessThan{Float64},MOI.EqualTo{Float64}}

@enum(_RowType, _ROW_TYPE_LESSTHAN, _ROW_TYPE_GREATERTHAN, _ROW_TYPE_EQUAL_TO,)

_row_type(::MOI.GreaterThan{Float64}) = _ROW_TYPE_GREATERTHAN
_row_type(::MOI.LessThan{Float64}) = _ROW_TYPE_LESSTHAN
_row_type(::MOI.EqualTo{Float64}) = _ROW_TYPE_EQUAL_TO

_bounds(s::MOI.GreaterThan{Float64}) = s.lower, Inf
_bounds(s::MOI.LessThan{Float64}) = -Inf, s.upper
_bounds(s::MOI.EqualTo{Float64}) = s.value, s.value

@enum(
    _BoundEnum,
    _BOUND_NONE,
    _BOUND_LESS_THAN,
    _BOUND_GREATER_THAN,
    _BOUND_LESS_AND_GREATER_THAN,
    _BOUND_EQUAL_TO,
)

@enum(_TypeEnum, _TYPE_CONTINUOUS, _TYPE_INTEGER, _TYPE_BINARY,)

"""
    _VariableInfo

A struct to store information about the variables.
"""
mutable struct _VariableInfo
    # We need to keep the index here because sometimes we call `LinearIndex`.
    index::MOI.VariableIndex
    # The variable name.
    name::String
    # The zero-indexed column in the cuOpt object.
    column::Int32
    # Storage to keep track of the variable bounds.
    bound::_BoundEnum
    lower::Float64
    upper::Float64
    # Track integrality
    type::_TypeEnum
    start::Union{Nothing,Float64}
    function _VariableInfo(
        index::MOI.VariableIndex,
        column::Int32,
        bound::_BoundEnum = _BOUND_NONE,
    )
        return new(
            index,
            "",
            column,
            bound,
            -Inf,
            Inf,
            _TYPE_CONTINUOUS,
            nothing,
        )
    end
end

function _variable_info_dict()
    return CleverDicts.CleverDict{MOI.VariableIndex,_VariableInfo}(
        x -> x.value,
        i -> MOI.VariableIndex(i),
    )
end

"""
    _ConstraintInfo

A struct to store information about the affine constraints.
"""
mutable struct _ConstraintInfo
    # The constraint name.
    name::String
    # The zero-indexed row in the cuOpt object.
    row::Int32
    # Storage to keep track of the constraint bounds.
    set::_RowType
    lower::Float64
    upper::Float64
end

function _ConstraintInfo(set::_SCALAR_SETS)
    lower, upper = _bounds(set)
    return _ConstraintInfo("", 0, _row_type(set), lower, upper)
end

struct _ConstraintKey
    value::Int64
end

function _constraint_info_dict()
    return CleverDicts.CleverDict{_ConstraintKey,_ConstraintInfo}(
        c -> c.value,
        i -> _ConstraintKey(i),
    )
end

"""
    _set(c::_ConstraintInfo)

Return the set associated with a constraint.
"""
function _set(c::_ConstraintInfo)
    if c.set == _ROW_TYPE_LESSTHAN
        return MOI.LessThan(c.upper)
    elseif c.set == _ROW_TYPE_GREATERTHAN
        return MOI.GreaterThan(c.lower)
    elseif c.set == _ROW_TYPE_INTERVAL
        return MOI.Interval(c.lower, c.upper)
    else
        @assert c.set == _ROW_TYPE_EQUAL_TO
        return MOI.EqualTo(c.lower)
    end
end

@enum(_OptimizeStatus, _OPTIMIZE_NOT_CALLED, _OPTIMIZE_OK, _OPTIMIZE_ERRORED)

mutable struct Optimizer <: MOI.AbstractOptimizer
    cuopt_problem::cuOptOptimizationProblem
    cuopt_settings::cuOptSolverSettings
    cuopt_solution::cuOptSolution

    name::String
    silent::Bool

    variable_info::typeof(_variable_info_dict())
    affine_constraint_info::typeof(_constraint_info_dict())

    raw_optimizer_attributes::Dict{String,Any}
    objective_sense::Union{Nothing,MOI.OptimizationSense}

    primal_solution::Vector{Float64}
    dual_solution::Vector{Float64}

    obj_value::Float64
    termination_status::Int32
    function Optimizer()
        model = new()
        model.cuopt_problem = cuOptOptimizationProblem(C_NULL)
        model.cuopt_settings = cuOptSolverSettings(C_NULL)
        model.cuopt_solution = cuOptSolution(C_NULL)

        model.name = ""
        model.silent = false

        model.variable_info = _variable_info_dict()
        model.affine_constraint_info = _constraint_info_dict()

        model.raw_optimizer_attributes = Dict{String,Any}(
            CUOPT_TIME_LIMIT => 3600,
            CUOPT_NUM_CPU_THREADS => 1,
            CUOPT_MIP_ABSOLUTE_GAP => 1e-10,
            CUOPT_MIP_RELATIVE_GAP => 1e-4,
        )

        model.objective_sense = nothing

        model.primal_solution = Float64[]
        model.dual_solution = Float64[]
        model.obj_value = 0.0
        model.termination_status = CUOPT_TERIMINATION_STATUS_NO_TERMINATION
        return model
    end
end

MOI.supports_incremental_interface(::Optimizer) = false

function MOI.is_empty(model::Optimizer)
    return length(model.variable_info) == 0 &&
           length(model.affine_constraint_info) == 0 &&
           model.objective_sense == nothing &&
           model.name == ""
end

"""
SolverName	Yes	No	No
SolverVersion	Yes	No	No
RawSolver	Yes	No	No
Name	Yes	Yes	Yes
Silent	Yes	Yes	Yes
TimeLimitSec	Yes	Yes	Yes
ObjectiveLimit	Yes	Yes	Yes
SolutionLimit	Yes	Yes	Yes
NodeLimit	Yes	Yes	Yes
RawOptimizerAttribute	Yes	Yes	Yes
NumberOfThreads	Yes	Yes	Yes
AbsoluteGapTolerance	Yes	Yes	Yes
RelativeGapTolerance	Yes	Yes	Yes
"""

function MOI.get(::Optimizer, ::MOI.SolverName)
    return "cuOpt"
end

function MOI.get(::Optimizer, ::MOI.SolverVersion)
    version_major = Ref{Int32}()
    version_minor = Ref{Int32}()
    version_patch = Ref{Int32}()
    ret = cuOptGetVersion(version_major, version_minor, version_patch)
    _check_ret(ret, "cuOptGetVersion")
    return "$(lpad(version_major[], 2, '0')).$(lpad(version_minor[], 2, '0')).$(version_patch[])"
end

#MOI.get(model::Optimizer, ::MOI.RawSolver) = model

function MOI.get(model::Optimizer, ::MOI.Name)
    return model.name
end

function MOI.get(model::Optimizer, ::MOI.Silent)
    return model.silent
end

function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    if param.name in keys(model.raw_optimizer_attributes)
        return model.raw_optimizer_attributes[param.name]
    end
    return nothing
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    return MOI.get(model, MOI.RawOptimizerAttribute(CUOPT_TIME_LIMIT))
end
#MOI.get(model::Optimizer, ::MOI.ObjectiveLimit) = model.objective_limit
#MOI.get(model::Optimizer, ::MOI.SolutionLimit) = model.solution_limit
#MOI.get(model::Optimizer, ::MOI.NodeLimit) = model.node_limit

function MOI.get(model::Optimizer, ::MOI.NumberOfThreads)
    return MOI.get(model, MOI.RawOptimizerAttribute(CUOPT_NUM_CPU_THREADS))
end

function MOI.get(model::Optimizer, ::MOI.AbsoluteGapTolerance)
    return MOI.get(model, MOI.RawOptimizerAttribute(CUOPT_MIP_ABSOLUTE_GAP))
end

function MOI.get(model::Optimizer, ::MOI.RelativeGapTolerance)
    return MOI.get(model, MOI.RawOptimizerAttribute(CUOPT_MIP_RELATIVE_GAP))
end

function MOI.get(model::Optimizer, ::MOI.NumberOfVariables)
    return length(model.variable_info)
end

const _TerminationStatusMap = Dict(
    CUOPT_TERIMINATION_STATUS_NO_TERMINATION =>
        (MOI.OPTIMIZE_NOT_CALLED, "cuOptModelStatusNotset"),
    CUOPT_TERIMINATION_STATUS_OPTIMAL =>
        (MOI.OPTIMAL, "cuOptModelStatusOptimal"),
    CUOPT_TERIMINATION_STATUS_INFEASIBLE =>
        (MOI.INFEASIBLE, "cuOptModelStatusInfeasible"),
    CUOPT_TERIMINATION_STATUS_UNBOUNDED =>
        (MOI.DUAL_INFEASIBLE, "cuOptModelStatusUnbounded"),
    CUOPT_TERIMINATION_STATUS_ITERATION_LIMIT =>
        (MOI.ITERATION_LIMIT, "cuOptModelStatusIterationLimit"),
    CUOPT_TERIMINATION_STATUS_TIME_LIMIT =>
        (MOI.TIME_LIMIT, "cuOptModelStatusTimeLimit"),
    CUOPT_TERIMINATION_STATUS_NUMERICAL_ERROR =>
        (MOI.NUMERICAL_ERROR, "cuOptModelStatusNumericalError"),
    CUOPT_TERIMINATION_STATUS_PRIMAL_FEASIBLE =>
        (MOI.OTHER_LIMIT, "cuOptModelStatusPrimalFeasible"),
    CUOPT_TERIMINATION_STATUS_CONCURRENT_LIMIT =>
        (MOI.OTHER_ERROR, "cuOptModelStatusConcurrent"),
    CUOPT_TERIMINATION_STATUS_FEASIBLE_FOUND =>
        (MOI.OTHER_LIMIT, "cuOptModelStatusFeasible"),
)

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    return _TerminationStatusMap[model.termination_status][1]
end

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    return _TerminationStatusMap[model.termination_status][2]
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveValue)
    return model.obj_value
end

function MOI.get(model::Optimizer, ::MOI.DualObjectiveValue)
    is_mip = Ref{Int32}()
    ret = cuOptIsMIP(model.cuopt_problem, is_mip)
    _check_ret(ret, "cuOptIsMIP")
    if is_mip[] == 1
        return NaN
    else
        dual_obj = Ref{Float64}()
        ret = cuOptGetDualObjectiveValue(model.cuopt_solution, dual_obj)
        _check_ret(ret, "cuOptGetDualObjectiveValue")
        return dual_obj[]
    end
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveBound)
    is_mip = Ref{Int32}()
    ret = cuOptIsMIP(model.cuopt_problem, is_mip)
    _check_ret(ret, "cuOptIsMIP")
    if is_mip[] == 1
        obj_bound = Ref{Float64}()
        ret = cuOptGetSolutionBound(model.cuopt_solution, obj_bound)
        _check_ret(ret, "cuOptGetSolutionBound")
        return obj_bound[]
    else
        # For now return the dual objective value for LP problems.
        dual_obj = Ref{Float64}()
        ret = cuOptGetDualObjectiveValue(model.cuopt_solution, dual_obj)
        _check_ret(ret, "cuOptGetDualObjectiveValue")
        return dual_obj[]
    end
end

function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    solve_time = Ref{Float64}()
    ret = cuOptGetSolveTime(model.cuopt_solution, solve_time)
    _check_ret(ret, "cuOptGetSolveTime")
    return solve_time[]
end

function MOI.get(model::Optimizer, ::MOI.ListOfVariableAttributesSet)
    ret = MOI.AbstractVariableAttribute[]
    found_name, found_start = false, false
    for info in values(model.variable_info)
        if !found_name && !isempty(info.name)
            push!(ret, MOI.VariableName())
            found_name = true
        end
        if !found_start && info.start !== nothing
            push!(ret, MOI.VariablePrimalStart())
            found_start = true
        end
        if found_start && found_name
            return ret
        end
    end
    return ret
end

function MOI.get(model::Optimizer, ::MOI.ListOfModelAttributesSet)
    if MOI.is_empty(model)
        return Any[]
    end
    attributes = Any[]
    if model.objective_sense !== nothing
        push!(attributes, MOI.ObjectiveSense())
    end
    # if model.is_objective_set
    #     F = MOI.get(model, MOI.ObjectiveFunctionType())
    #     push!(attributes, MOI.ObjectiveFunction{F}())
    # end
    if MOI.get(model, MOI.Name()) != ""
        push!(attributes, MOI.Name())
    end

    return attributes
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintAttributesSet{F,S},
) where {S,F}
    if F == MOI.VariableIndex
        # Does not support ConstraintName
        return MOI.AbstractConstraintAttribute[]
    end
    for index in MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
        if !isempty(MOI.get(model, MOI.ConstraintName(), index))
            return MOI.AbstractConstraintAttribute[MOI.ConstraintName()]
        end
    end
    return MOI.AbstractConstraintAttribute[]
end

function _check_ret(ret::Int32, msg::String)
    if ret == CUOPT_SUCCESS
        return
    elseif ret == CUOPT_INVALID_ARGUMENT
        error("Invalid argument in $msg.")
    elseif ret == CUOPT_MPS_FILE_ERROR
        error("MPS file open error in $msg.")
    elseif ret == CUOPT_MPS_PARSE_ERROR
        error("MPS parse error in $msg.")
    elseif ret == CUOPT_VALIDATION_ERROR
        error("Validation error in $msg.")
    elseif ret == CUOPT_OUT_OF_MEMORY
        error("Out of memory in $msg.")
    elseif ret == CUOPT_RUNTIME_ERROR
        error("Runtime error in $msg.")
    end
end

function MOI.get(
    model::Optimizer,
    attr::MOI.VariablePrimal,
    x::MOI.VariableIndex,
)
    #MOI.check_result_index_bounds(model, attr)
    return model.primal_solution[column(model, x)+1]
end

function MOI.get(model::Optimizer, ::MOI.ResultCount)
    status = model.termination_status
    if status == CUOPT_TERIMINATION_STATUS_OPTIMAL ||
       status == CUOPT_TERIMINATION_STATUS_PRIMAL_FEASIBLE ||
       status == CUOPT_TERIMINATION_STATUS_FEASIBLE_FOUND ||
       status == CUOPT_TERIMINATION_STATUS_UNBOUNDED
        return 1
    end
    return 0
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    status = model.termination_status
    if attr.result_index != 1
        return MOI.NO_SOLUTION
    elseif status == CUOPT_TERIMINATION_STATUS_OPTIMAL ||
           status == CUOPT_TERIMINATION_STATUS_PRIMAL_FEASIBLE ||
           status == CUOPT_TERIMINATION_STATUS_FEASIBLE_FOUND
        return MOI.FEASIBLE_POINT
    elseif status == CUOPT_TERIMINATION_STATUS_INFEASIBLE
        return MOI.INFEASIBLE_POINT
    end
    return MOI.NO_SOLUTION
end

function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    status = model.termination_status
    if attr.result_index != 1
        return MOI.NO_SOLUTION
    elseif status == CUOPT_TERIMINATION_STATUS_OPTIMAL
        return MOI.FEASIBLE_POINT
    end
    return MOI.NO_SOLUTION
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveSense)
    if model.objective_sense !== nothing
        return model.objective_sense
    end
    return MOI.MIN_SENSE
end

function MOI.set(model::Optimizer, ::MOI.Name, name::String)
    return model.name = name
end

function MOI.set(model::Optimizer, ::MOI.Silent, silent::Bool)
    return model.silent = silent
end

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    return model.raw_optimizer_attributes[param.name] = value
end

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, time_limit::Real)
    return MOI.set(
        model,
        MOI.RawOptimizerAttribute(CUOPT_TIME_LIMIT),
        time_limit,
    )
end

#MOI.set(model::Optimizer, ::MOI.ObjectiveLimit, objective_limit::Float64) = model.objective_limit = objective_limit
#MOI.set(model::Optimizer, ::MOI.SolutionLimit, solution_limit::Int) = model.solution_limit = solution_limit
#MOI.set(model::Optimizer, ::MOI.NodeLimit, node_limit::Int) = model.node_limit = node_limit

function MOI.set(
    model::Optimizer,
    ::MOI.NumberOfThreads,
    number_of_threads::Int,
)
    return MOI.set(
        model,
        MOI.RawOptimizerAttribute(CUOPT_NUM_CPU_THREADS),
        number_of_threads,
    )
end

function MOI.set(
    model::Optimizer,
    ::MOI.AbsoluteGapTolerance,
    absolute_gap_tolerance::Float64,
)
    return MOI.set(
        model,
        MOI.RawOptimizerAttribute(CUOPT_MIP_ABSOLUTE_GAP),
        absolute_gap_tolerance,
    )
end

function MOI.set(
    model::Optimizer,
    ::MOI.RelativeGapTolerance,
    relative_gap_tolerance::Float64,
)
    return MOI.set(
        model,
        MOI.RawOptimizerAttribute(CUOPT_MIP_RELATIVE_GAP),
        relative_gap_tolerance,
    )
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveSense,
    sense::MOI.OptimizationSense,
)
    if !(sense in (MOI.MIN_SENSE, MOI.MAX_SENSE))
        throw(MOI.UnsupportedAttribute(MOI.ObjectiveSense()))
    end

    return model.objective_sense = sense
end

# features that cuOpt supports
MOI.supports(::Optimizer, ::MOI.Name) = true
MOI.supports(::Optimizer, ::MOI.RawSolver) = false
MOI.supports(::Optimizer, ::MOI.Silent) = true
MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveLimit) = false
MOI.supports(::Optimizer, ::MOI.SolutionLimit) = false
MOI.supports(::Optimizer, ::MOI.NodeLimit) = false
MOI.supports(::Optimizer, ::MOI.RawOptimizerAttribute) = true
MOI.supports(::Optimizer, ::MOI.NumberOfThreads) = true
MOI.supports(::Optimizer, ::MOI.AbsoluteGapTolerance) = true
MOI.supports(::Optimizer, ::MOI.RelativeGapTolerance) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
MOI.supports(::Optimizer, ::MOI.ConstraintFunction) = true
MOI.supports(::Optimizer, ::MOI.DualObjectiveValue) = false

function MOI.empty!(model::Optimizer)
    model.name = ""
    # attribute silent should not be cleared
    model.variable_info = _variable_info_dict()
    model.affine_constraint_info = _constraint_info_dict()
    model.objective_sense = nothing
    model.primal_solution = Float64[]
    model.dual_solution = Float64[]
    model.obj_value = NaN
    model.termination_status = CUOPT_TERIMINATION_STATUS_NO_TERMINATION

    if model.cuopt_problem != C_NULL
        a = Ref(model.cuopt_problem)
        cuOptDestroyProblem(a)
        model.cuopt_problem = C_NULL
    end

    if model.cuopt_settings != C_NULL
        a = Ref(model.cuopt_settings)
        cuOptDestroySolverSettings(a)
        model.cuopt_settings = C_NULL
    end

    if model.cuopt_solution != C_NULL
        a = Ref(model.cuopt_solution)
        cuOptDestroySolution(a)
        model.cuopt_solution = C_NULL
    end

    return
end

# supports constraints on variables
function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{<:_SCALAR_SETS},
)
    return true
end

# supports integrality constraints on variables
function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{<:Union{MOI.ZeroOne,MOI.Integer}},
)
    return true
end

# does not support semicontinuous and semiinteger constraints on variables
function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{S},
) where {S<:Union{MOI.Semicontinuous{Float64},MOI.Semiinteger}}
    return false
end

# supports constraints on affine functions
function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarAffineFunction{Float64}},
    ::Type{<:_SCALAR_SETS},
)
    return true
end

function MOI.supports(
    ::Optimizer,
    ::Union{MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}},
)
    return true
end

function _check_input_data(dest::Optimizer, src::MOI.ModelLike)
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        if !MOI.supports_constraint(dest, F, S)
            throw(
                MOI.UnsupportedConstraint{F,S}(
                    "cuOpt does not support constraints of type $F-in-$S.",
                ),
            )
        end
        for attr in MOI.get(src, MOI.ListOfConstraintAttributesSet{F,S}())
            if attr in (
                MOI.ConstraintName(),
                MOI.ConstraintPrimalStart(),
                MOI.ConstraintDualStart(),
            )
                continue
            else
                throw(MOI.UnsupportedAttribute(attr))
            end
        end
    end
    for attr in MOI.get(src, MOI.ListOfModelAttributesSet())
        if attr in (
            MOI.Name(),
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            MOI.ObjectiveSense(),
        )
            continue
        else
            throw(MOI.UnsupportedAttribute(attr))
        end
    end
    for attr in MOI.get(src, MOI.ListOfVariableAttributesSet())
        if attr in (MOI.VariableName(), MOI.VariablePrimalStart())
            continue
        else
            throw(MOI.UnsupportedAttribute(attr))
        end
    end
    return
end

function _info(model::Optimizer, key::MOI.VariableIndex)
    info = get(model.variable_info, key, nothing)
    if info === nothing
        return throw(MOI.InvalidIndex(key))
    end
    return info
end

"""
    column(model::Optimizer, x::MOI.VariableIndex)

Return the 0-indexed column associated with `x` in `model`.
"""
column(model::Optimizer, x::MOI.VariableIndex) = _info(model, x).column

function _info(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},
)
    var_index = MOI.VariableIndex(c.value)
    info = get(model.variable_info, var_index, nothing)
    if info === nothing
        return throw(MOI.InvalidIndex(c))
    end
    return info
end

"""
    column(
        model::Optimizer,
        c::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},
    )

Return the 0-indexed column associated with the variable bounds `c` in `model`.
"""
function column(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},
)
    return _info(model, c).column
end

function _numcols(model::Optimizer)
    return Int32(length(model.variable_info))
end

function _numrows(model::Optimizer)
    return Int32(length(model.affine_constraint_info))
end

function _copy_to_columns(dest::Optimizer, src::MOI.ModelLike, mapping)
    x_src = MOI.get(src, MOI.ListOfVariableIndices())
    numcols = Int32(length(x_src))
    for (i, x) in enumerate(x_src)
        index = CleverDicts.add_item(
            dest.variable_info,
            _VariableInfo(MOI.VariableIndex(0), Int32(0)),
        )
        info = _info(dest, index)
        if MOI.supports(src, MOI.VariableName(), MOI.VariableIndex)
            info.name = MOI.get(src, MOI.VariableName(), x)
        end
        info.index = index
        info.column = Int32(i - 1)
        mapping[x] = index
    end
    return numcols
end

function _add_bounds(::Vector{Float64}, ub, i, s::MOI.LessThan{Float64})
    ub[i] = s.upper
    return
end

function _add_bounds(lb, ::Vector{Float64}, i, s::MOI.GreaterThan{Float64})
    lb[i] = s.lower
    return
end

function _add_bounds(lb, ub, i, s::MOI.EqualTo{Float64})
    lb[i], ub[i] = s.value, s.value
    return
end

function _throw_if_existing_lower(
    info::_VariableInfo,
    ::S,
) where {S<:MOI.AbstractSet}
    if info.bound == _BOUND_LESS_AND_GREATER_THAN
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Float64},S}(info.index))
    elseif info.bound == _BOUND_GREATER_THAN
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Float64},S}(info.index))
    elseif info.bound == _BOUND_EQUAL_TO
        throw(MOI.LowerBoundAlreadySet{MOI.EqualTo{Float64},S}(info.index))
    end
    return
end

function _throw_if_existing_upper(
    info::_VariableInfo,
    ::S,
) where {S<:MOI.AbstractSet}
    if info.bound == _BOUND_LESS_AND_GREATER_THAN
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Float64},S}(info.index))
    elseif info.bound == _BOUND_LESS_THAN
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Float64},S}(info.index))
    elseif info.bound == _BOUND_EQUAL_TO
        throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{Float64},S}(info.index))
    end
    return
end

function _update_info(info::_VariableInfo, s::MOI.GreaterThan{Float64})
    _throw_if_existing_lower(info, s)
    if info.bound == _BOUND_LESS_THAN
        info.bound = _BOUND_LESS_AND_GREATER_THAN
    else
        info.bound = _BOUND_GREATER_THAN
    end
    info.lower = s.lower
    return
end

function _update_info(info::_VariableInfo, s::MOI.LessThan{Float64})
    _throw_if_existing_upper(info, s)
    if info.bound == _BOUND_GREATER_THAN
        info.bound = _BOUND_LESS_AND_GREATER_THAN
    else
        info.bound = _BOUND_LESS_THAN
    end
    info.upper = s.upper
    return
end

function _update_info(info::_VariableInfo, s::MOI.EqualTo{Float64})
    _throw_if_existing_lower(info, s)
    _throw_if_existing_upper(info, s)
    info.bound = _BOUND_EQUAL_TO
    info.lower = s.value
    info.upper = s.value
    return
end

_constraints(src, F, S) = MOI.get(src, MOI.ListOfConstraintIndices{F,S}())

function _extract_bound_data(
    dest::Optimizer,
    src::MOI.ModelLike,
    mapping,
    collower::Vector{Float64},
    colupper::Vector{Float64},
    ::Type{S},
) where {S}
    for c_index in _constraints(src, MOI.VariableIndex, S)
        f = MOI.get(src, MOI.ConstraintFunction(), c_index)
        s = MOI.get(src, MOI.ConstraintSet(), c_index)
        new_f = mapping[f]
        info = _info(dest, new_f)
        _add_bounds(collower, colupper, info.column + 1, s)
        _update_info(info, s)
        mapping[c_index] = MOI.ConstraintIndex{MOI.VariableIndex,S}(new_f.value)
    end
    return
end

_add_sizehint!(vec, n) = sizehint!(vec, length(vec) + n)

function _extract_row_data(
    dest::Optimizer,
    src::MOI.ModelLike,
    mapping,
    rowlower::Vector{Float64},
    rowupper::Vector{Float64},
    constraint_matrix_row_offsets::Vector{Int32},
    constraint_matrix_column_indices::Vector{Int32},
    constraint_matrix_coefficients::Vector{Float64},
    ::Type{S},
) where {S}
    F = MOI.ScalarAffineFunction{Float64}
    row = length(rowlower) + 1
    list = _constraints(src, MOI.ScalarAffineFunction{Float64}, S)
    numrows = length(list)
    _add_sizehint!(rowlower, numrows)
    _add_sizehint!(rowupper, numrows)
    n_terms = 0
    fs = Array{F}(undef, numrows)
    for (i, c_index) in enumerate(list)
        f = MOI.get(src, MOI.ConstraintFunction(), c_index)
        fs[i] = f
        set = MOI.get(src, MOI.ConstraintSet(), c_index)
        l, u = _bounds(set)
        push!(rowlower, l - f.constant)
        push!(rowupper, u - f.constant)
        n_terms += length(f.terms)
        key = CleverDicts.add_item(
            dest.affine_constraint_info,
            _ConstraintInfo(set),
        )
        info = dest.affine_constraint_info[key]
        info.row = Int32(length(dest.affine_constraint_info) - 1)
        if MOI.supports(src, MOI.ConstraintName(), MOI.ConstraintIndex{F,S})
            info.name = MOI.get(src, MOI.ConstraintName(), c_index)
        end
        mapping[c_index] = MOI.ConstraintIndex{F,S}(key.value)
    end

    nnz = constraint_matrix_row_offsets[end]
    for f in fs
        for term in f.terms
            push!(
                constraint_matrix_column_indices,
                Int32(mapping[term.variable].value) - 1,
            )
            push!(constraint_matrix_coefficients, term.coefficient)
            nnz += 1
        end
        row += 1
        push!(constraint_matrix_row_offsets, nnz)
    end
    return
end

function _get_objective_data(
    dest::Optimizer,
    src::MOI.ModelLike,
    mapping,
    numcol::Int32,
)
    sense = MOI.get(src, MOI.ObjectiveSense())
    if !(sense in (MOI.MIN_SENSE, MOI.MAX_SENSE))
        throw(MOI.UnsupportedAttribute(MOI.ObjectiveSense()))
    end

    objective_sense = sense == MOI.MIN_SENSE ? CUOPT_MINIMIZE : CUOPT_MAXIMIZE

    objective_coefficients = zeros(Float64, numcol)
    F = MOI.get(src, MOI.ObjectiveFunctionType())
    f_obj = MOI.get(src, MOI.ObjectiveFunction{F}())

    for term in f_obj.terms
        objective_coefficients[mapping[term.variable].value] += term.coefficient
    end

    objective_offset = f_obj.constant

    return objective_sense, objective_offset, objective_coefficients
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    MOI.empty!(dest)
    _check_input_data(dest, src)

    mapping = MOI.Utilities.IndexMap()
    numcol = _copy_to_columns(dest, src, mapping)

    collower, colupper = fill(-Inf, numcol), fill(Inf, numcol)
    rowlower, rowupper = Float64[], Float64[]

    constraint_matrix_row_offsets = Int32[]
    constraint_matrix_column_indices = Int32[]
    constraint_matrix_coefficients = Float64[]

    push!(constraint_matrix_row_offsets, 0)
    for S in
        (MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.EqualTo{Float64})
        _extract_bound_data(dest, src, mapping, collower, colupper, S)
        _extract_row_data(
            dest,
            src,
            mapping,
            rowlower,
            rowupper,
            constraint_matrix_row_offsets,
            constraint_matrix_column_indices,
            constraint_matrix_coefficients,
            S,
        )
    end

    numrow = Int32(length(rowlower))
    # Extract integrality constraints
    var_type = fill(Cchar('C'), numcol)
    has_integrality = false
    for ci in _constraints(src, MOI.VariableIndex, MOI.ZeroOne)
        info = _info(dest, ci)
        info.type = _TYPE_BINARY
        var_type[info.column+1] = Cchar('I')
        new_x = mapping[MOI.VariableIndex(ci.value)]
        mapping[ci] = typeof(ci)(new_x.value)
        collower[info.column+1] = ceil(max(0.0, collower[info.column+1]))
        colupper[info.column+1] = floor(min(1.0, colupper[info.column+1]))
        has_integrality = true
    end
    for ci in _constraints(src, MOI.VariableIndex, MOI.Integer)
        info = _info(dest, ci)
        info.type = _TYPE_INTEGER
        var_type[info.column+1] = Cchar('I')
        new_x = mapping[MOI.VariableIndex(ci.value)]
        mapping[ci] = typeof(ci)(new_x.value)
        has_integrality = true
    end

    objective_sense, objective_offset, objective_coefficients =
        _get_objective_data(dest, src, mapping, numcol)

    ref_problem = Ref{cuOptOptimizationProblem}()
    ret = cuOptCreateRangedProblem(
        numrow,
        numcol,
        objective_sense,
        objective_offset,
        objective_coefficients,
        constraint_matrix_row_offsets,
        constraint_matrix_column_indices,
        constraint_matrix_coefficients,
        rowlower,
        rowupper,
        collower,
        colupper,
        var_type,
        ref_problem,
    )
    _check_ret(ret, "cuOptCreateRangedProblem")
    dest.cuopt_problem = ref_problem[]

    ref_settings = Ref{cuOptSolverSettings}()
    ret = cuOptCreateSolverSettings(ref_settings)
    _check_ret(ret, "cuOptCreateSolverSettings")
    dest.cuopt_settings = ref_settings[]

    # Set all raw optimizer attributes
    for (name, value) in dest.raw_optimizer_attributes
        ret = cuOptSetParameter(dest.cuopt_settings, name, string(value))
        _check_ret(ret, "cuOptSetParameter($name, $value)")
    end

    # Override log info if silent is set
    if dest.silent
        ret =
            cuOptSetParameter(dest.cuopt_settings, "log_to_console", string(0))
        _check_ret(ret, "cuOptSetParameter(log_to_console, 0)")
        ret = cuOptSetParameter(dest.cuopt_settings, "log_file", "")
        _check_ret(ret, "cuOptSetParameter(log_file, )")
    end

    return mapping
end

function MOI.optimize!(model::Optimizer)
    ref_solution = Ref{cuOptSolution}()
    ret = cuOptSolve(model.cuopt_problem, model.cuopt_settings, ref_solution)
    _check_ret(ret, "cuOptSolve")
    model.cuopt_solution = ref_solution[]

    status = Ref{Int32}()
    ret = cuOptGetTerminationStatus(model.cuopt_solution, status)
    _check_ret(ret, "cuOptGetTerminationStatus")
    model.termination_status = status[]

    obj_value = zeros(1)
    ret = cuOptGetObjectiveValue(model.cuopt_solution, obj_value)
    _check_ret(ret, "cuOptGetObjectiveValue")
    model.obj_value = obj_value[1]

    model.primal_solution = zeros(_numcols(model))
    model.dual_solution = zeros(_numrows(model))

    ret = cuOptGetPrimalSolution(model.cuopt_solution, model.primal_solution)
    _check_ret(ret, "cuOptGetPrimalSolution")

    is_mip = Ref{Int32}()
    ret = cuOptIsMIP(model.cuopt_problem, is_mip)
    _check_ret(ret, "cuOptIsMIP")
    if is_mip[] == 0
        ret = cuOptGetDualSolution(model.cuopt_solution, model.dual_solution)
        _check_ret(ret, "cuOptGetDualSolution")
    end

    return
end
