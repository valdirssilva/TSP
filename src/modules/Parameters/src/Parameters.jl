module Parameters

using JuMP
using DelimitedFiles

Base.@kwdef mutable struct ExperimentParameters
    approach::String = MIP_solver # {MIP_solver}
    solver::String = Gurobi # {Gurobi, CPLEX}
    total_time_limit::Float64 = 1e6 # Total time limit for approach in seconds
    MIP_gap_tolerance::Float64 = 1e8
    integer_feasibility_tolerance::Float64 = 1e6
    number_of_threads::Int64 = 1
    screen_output::Int64 = 1
    approach_HC::String = random
    approach_LS::String = bestImprovementReinsertion
    approach_improvemente::String = Reinsertion
    rf_max_time::Float64 = 1800.0 # Time limit for relax-and-fix in seconds
    rf_restricted_model_max_time::Float64 = 10.0 # Time limit for RF restricted models in seconds
    rf_Number_of_Int::Float64 = 4 # Portion of number of periods in instance. rf_window_size = max{2; floor(number_of_periods / rf_window_size_factor)}
    rf_Number_of_fix::Float64 = 2  # Portion of number of window size. rf_step_size = ceil(window_size / rf_step_size_factor)
    
end

export ExperimentParameters, read_parameters,set_MIP_solver_parameters

function read_parameters(parameters_file::String)

    param_data = readdlm(parameters_file)

    approach = param_data[2,2]
    solver = param_data[3,2]
    total_time_limit = param_data[4,2]
    MIP_gap_tolerance = param_data[5,2]
    integer_feasibility_tolerance = param_data[6,2]
    number_of_threads = param_data[7,2]
    screen_output = param_data[8,2]
    approach_HC = param_data[9,2]
    approach_LS = param_data[10,2]
    approach_improvemente = param_data[11,2]
    rf_max_time = param_data[13,2]
    rf_restricted_model_max_time = param_data[14,2]
    rf_Number_of_Int = param_data[15,2]
    rf_Number_of_fix = param_data[16,2]

    parameters = ExperimentParameters(
        approach = approach,
        solver = solver,
        total_time_limit = total_time_limit,
        MIP_gap_tolerance = MIP_gap_tolerance,
        integer_feasibility_tolerance = integer_feasibility_tolerance,
        number_of_threads = number_of_threads,
        screen_output = screen_output,
        approach_HC = approach_HC,
        approach_LS = approach_LS,
        approach_improvemente = approach_improvemente,
        rf_max_time = rf_max_time,
        rf_restricted_model_max_time = rf_restricted_model_max_time,
        rf_Number_of_Int = rf_Number_of_Int,
        rf_Number_of_fix = rf_Number_of_fix
       
    )
    return parameters
end
function set_MIP_solver_parameters(model::Model, params::ExperimentParameters)
    
    if params.solver == "Gurobi"
        set_optimizer_attributes(model,
                                "TimeLimit" => params.total_time_limit,
                                "MIPGap" => params.MIP_gap_tolerance,
                                "IntFeasTol" => params.integer_feasibility_tolerance,
                                "Threads" => params.number_of_threads,
                                "LogToConsole" => params.screen_output)
    elseif params.solver == "CPLEX"
        set_optimizer_attributes(model,
                                "CPXPARAM_TimeLimit" => params.total_time_limit,
                                "CPXPARAM_MIP_Tolerances_MIPGap" => params.MIP_gap_tolerance,
                                "CPXPARAM_MIP_Tolerances_Integrality" => params.integer_feasibility_tolerance,
                                "CPX_PARAM_THREADS" => params.number_of_threads,
                                "CPXPARAM_ScreenOutput" => params.screen_output)
    end
end
function set_rf_parameters(model::Model, params::ExperimentParameters)
    if params.solver == "Gurobi"
        set_optimizer_attributes(model,
                                "TimeLimit" => params.rf_restricted_model_max_time,
                                "MIPGap" => params.MIP_gap_tolerance,
                                "IntFeasTol" => params.integer_feasibility_tolerance,
                                "Threads" => params.number_of_threads,
                                "LogToConsole" => params.screen_output)
    elseif params.solver == "CPLEX"
        set_optimizer_attributes(model,
                                "CPXPARAM_TimeLimit" => params.rf_restricted_model_max_time,
                                "CPXPARAM_MIP_Tolerances_MIPGap" => params.MIP_gap_tolerance,
                                "CPXPARAM_MIP_Tolerances_Integrality" => params.integer_feasibility_tolerance,
                                "CPX_PARAM_THREADS" => params.number_of_threads,
                                "CPXPARAM_ScreenOutput" => params.screen_output)
    end
end

function set_fo_parameters(model::Model, params::ExperimentParameters)
    if params.solver == "Gurobi"
        set_optimizer_attributes(model,
                                "TimeLimit" => params.fo_restricted_model_max_time,
                                "MIPGap" => params.MIP_gap_tolerance,
                                "IntFeasTol" => params.integer_feasibility_tolerance,
                                "Threads" => params.number_of_threads,
                                "LogToConsole" => params.screen_output)
    elseif params.solver == "CPLEX"
        set_optimizer_attributes(model,
                                "CPXPARAM_TimeLimit" => params.fo_restricted_model_max_time,
                                "CPXPARAM_MIP_Tolerances_MIPGap" => params.MIP_gap_tolerance,
                                "CPXPARAM_MIP_Tolerances_Integrality" => params.integer_feasibility_tolerance,
                                "CPX_PARAM_THREADS" => params.number_of_threads,
                                "CPXPARAM_ScreenOutput" => params.screen_output)
    end
end



end # module