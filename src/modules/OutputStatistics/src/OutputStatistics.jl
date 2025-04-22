module OutputStatistics

using Data
using Parameters
using Dates

Base.@kwdef mutable struct StatisticsData
    date::String = string(Dates.today())
    approach::String = ""
    best_LB::Float64 = -1e8
    best_UB::Float64 = 1e8
    gap::Float64 = 100.0
    total_time::Float64 = 0.0
    node_count::Int64 = 0
    sol_status = 0
    UB_HC::Float64 = -1e8
    UB_LS::Float64 = -1e8
    LS_time::Float64 = 0.0
    HC_time::Float64 = 0.0
    improvement::Float64 = 0.0
    rf_UB::Float64 = 1e8
    rf_iterations::Int64 = 0
    rf_time::Float64 = 0.0
    rf_backtrack_attempts::Int64 = 0
    rf_backtrack__rounds::Int64 = 0
    fo_UB::Float64 = 1e8
    fo_improvement::Float64 = 0.0
    fo_iterations::Int64 = 0
    fo_rounds::Int64 = 0
    fo_time::Float64 = 0.0

    
end

mutable struct StdFormModelSolution
    primal_bound::Float64
    dual_bound::Float64
    x::Array{Float64}
    route::Array{Int64}
    status
end



export StatisticsData, StdFormModelSolution, setup_stats_file, print_stats, print_tur

function init_std_form_solution(data_type::InstanceTypeData)
    
    N = data_type.DIMENSION
    
    primal_bound = 1e8
    dual_bound = -1e8
    x = Array{Float64}(undef,N,N)
    route = Array{Int64}(undef,N)
    fill!(x, 0.0)

    status = 0

    solution = StdFormModelSolution(primal_bound,
                                    dual_bound,
                                    x,
                                    route,
                                    status)
            
    return solution

end

function setup_stats_file(params::ExperimentParameters, inputlist_file::String, parameters_file::String)

    inputlist_file = split(inputlist_file, "/")
    inputclass = "_" * inputlist_file[end]
    params_file = split(parameters_file, "/")
    params_file = "_" * params_file[end]    
    output_file_path = "src/outputFiles/tables/$(params.approach)/"
    if params.approach == "MIP_solver"   || params.approach == "RF" 
        output_file_path = output_file_path 

        output_file_path = output_file_path * "/$(params.solver)/$(Int64(params.total_time_limit))s/"
    elseif params.approach == "HC" 
        output_file_path = output_file_path 

        output_file_path = output_file_path * "/$(params.approach_HC)/"
    elseif params.approach == "HC_LS"

        output_file_path = output_file_path 

        if params.approach_LS == "bestImprovement" 
           
            output_file_path = output_file_path * "/$(params.approach_HC)/$(params.approach_improvemente)/$(params.approach_LS)/"
        else
            output_file_path = output_file_path * "/$(params.approach_HC)/$(params.approach_LS)/"
        end

    elseif params.approach == "ILS" || params.approach == "MultiStart" || params.approach == "SA" || params.approach == "VNS"|| params.approach == "GRASP" || params.approach == "TABU"

        output_file_path = output_file_path 

        output_file_path = output_file_path * "/$(params.approach)/"

    end

    mkpath(output_file_path)

    date_time = Dates.now()
    time_stamp = string(Dates.today(), "-",
        Dates.hour(date_time), "h",
        Dates.minute(date_time), "m",
        Dates.second(date_time), "s")
    output_file = output_file_path * time_stamp #= * inputclass  =#*params_file

    out = open(output_file,"w")

    println(out, "Statistics for TSP ")
    println(out, "Date: ", string(Dates.today()))

    if params.approach == "MIP_solver"
        println(out, "Approach: Formulation with MIP solver")
        println(out, "Solver: ", params.solver)
        println(out, "Time limit: ", params.total_time_limit)
        println(out, "MIP gap tolerance: ", params.MIP_gap_tolerance)
        println(out, "Integer feasibility tolerance: ", params.integer_feasibility_tolerance)
        println(out, "Number of threads: ", params.number_of_threads)
        println(out, "Screen output: ", params.screen_output)

        print(out, "\nInstance & N & LB & UB & gap & status & total_time \\\\")
    elseif params.approach == "HC" 
        println(out, "Approach: heuristic", params.approach)
        println(out, "HC Approach: ", params.approach_HC)

        print(out, "\nInstance & N & UB & total_time \\\\")
    elseif params.approach == "HC_LS"
        println(out, "HC Approach: ", params.approach)
        println(out, "HC Approach: ", params.approach_HC)
        println(out, "LS approach: ", params.approach_LS)
        print(out, "\nInstance & N & UB_hC & HC_time & UB_LS & LS_time & total_time \\\\")  
    elseif params.approach == "ILS" || params.approach == "MultiStart" || params.approach == "SA" || params.approach == "VNS"|| params.approach == "GRASP"|| params.approach == "TABU"
        println(out, "Approach: ", params.approach)
        print(out, "\nInstance & N & UB_hC & HC_time & UB_ILS & ILS_time & total_time \\\\") 
    end

    close(out)

    return output_file
end
function print_stats( data_type::InstanceTypeData,  params::ExperimentParameters, stats::StatisticsData, stats_file::String)
    out = open(stats_file,"a")

    print(out, "\n", data_type.instName, " & ",
    data_type.DIMENSION, " & ")

    if params.approach == "MIP_solver" || params.approach == "RF"
        print(out, round(stats.best_LB, digits = 2), " & ",
        round(stats.best_UB, digits = 2), " & ",
        round(stats.gap, digits = 4), " & ",
        stats.sol_status, " & ",
        round(stats.total_time, digits = 2), " \\\\")
    elseif params.approach == "HC" 
        print(out,stats.UB_HC, " & ",
        round(stats.HC_time, digits = 6), " \\\\") 
    elseif params.approach == "HC_LS" || params.approach == "ILS" || params.approach == "MultiStart" || params.approach == "SA"|| params.approach == "VNS" || params.approach == "GRASP"|| params.approach == "TABU"
        print(out,stats.UB_HC, " & ",
        round(stats.HC_time, digits = 6), " & ", 
        stats.UB_LS, " & ",
        round(stats.improvement, digits = 6), " & ",
        round(stats.LS_time, digits = 6), " & ",
        round(stats.total_time, digits = 6), " \\\\")      
    end

    close(out)

    return
end
function print_tur( data::InstanceData, data_type::InstanceTypeData, solution::StdFormModelSolution)




    N = 1:data_type.DIMENSION


  


    aux = 1
    in = 1
    while aux <= data_type.DIMENSION + 1
        i = in
        for j in N
            if solution.x[i, j] #= * data.c[i,j] =# > 0.0001
                print("$i ->  ")
                in = j
                break
            end
        end

        aux += 1
    end
   

    return
end
end # module