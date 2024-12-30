push!(LOAD_PATH, "src/modules/")

__precompile__()

using Pkg
Pkg.activate(".")

using JuMP
using Gurobi
using CPUTime
#using CPLEX
using DelimitedFiles
using Combinatorics
using Random



import Data
import OutputStatistics
import Parameters
import Formulation
import Constructive
import BestImprovement
import Disturbance
import MetaHeuristics

# >>>>> INSTRUCTIONS <<<<<
# Args:                        [1]              [2]             [3]                           
# To run the code: julia clsp_sc.jl <inputListFile> <parametersFile>
#
# EXAMPLE
# julia src/TSP.jl src/inputFiles/inputList_test.txt src/paramFiles/parameters
function main(ARGS)


    # Read inputlist file
    inputlist_file = String(ARGS[1])
    input = readdlm(inputlist_file)

    # Read parameters file
    parameters_file = String(ARGS[2])
    params = Parameters.read_parameters(parameters_file)

   

    # Set Gurobi environment
    if params.solver == "Gurobi"
        GRB_ENV = Gurobi.Env()
    else
        GRB_ENV = 0
    end

    # Setup statistics output file
    output_file = OutputStatistics.setup_stats_file(params, inputlist_file, parameters_file)

    # Get number of instances to run
    num_inst = input[1]

    # Run all num_inst instances in the input list file
    for inst in 1:num_inst

        # Initialize statistics data structure
        stats = OutputStatistics.StatisticsData()

        stats.approach = params.approach

        # Get the name of the file containing the instance data
        instance_file = String(input[inst + 1])   
        
        # Read instance data

        data_type = Data.readTypeData(instance_file)
        #data = Data.read_instance(instance_file)

        if data_type.COORD_Type == "EXPLICIT"
            #println("Weights are listed explicitly in the corresponding section")
            data = Data.readDataEXPLICIT(instance_file, data_type)
        else 
 
            data = Data.readDataXY(instance_file, data_type)
        end

        
        println("\nSTART >>>>>>>>>>>>>>>>>>>>>>>>>>")
        println("\nRunning instance $inst: ", data_type.instName)

        # Create an empty model
        model = Model()
        
        if params.approach == "MIP_solver"

            if params.solver == "Gurobi"
                model = Model(() -> Gurobi.Optimizer(GRB_ENV))
            elseif params.solver == "CPLEX"
                model = Model(() -> CPLEX.Optimizer())
            end

            start_time = time_ns()
            # Setup std formulation solution data structure
            solution = OutputStatistics.init_std_form_solution(data_type)

            # Formulate the model
            Formulation.create_stdform_model!(data,data_type, model)
            #
            Parameters.set_MIP_solver_parameters(model, params)

            #solver
            Formulation.solve_stdform_model!(model,data, data_type, solution, stats)

            finish_time = time_ns()
            total_time = (finish_time - start_time) * 1e-9
            stats.total_time = total_time

            #print_tur in screen

            println("tour", solution.route)
            OutputStatistics.print_tur(data,data_type, solution)

            println("\nFINISH <<<<<<<<<<<<<<<<<<<<<<<<<")
        elseif params.approach == "HC" || params.approach == "HC_LS"
            if params.approach_HC == "VMP"

                start_time = time_ns()

                # Setup std formulation solution data structure
                solution = OutputStatistics.init_std_form_solution(data_type)

                #solver for VMP

                Constructive.solve_VMP!(data, data_type, solution, stats)

                finish_HC_time = time_ns()
                total_time_hc = (finish_HC_time - start_time) * 1e-9
                stats.HC_time = total_time_hc


            elseif params.approach_HC == "BN"

                start_time = time_ns()

                # Setup std formulation solution data structure
                solution = OutputStatistics.init_std_form_solution(data_type)

                #solver for VMP

                Constructive.solve_BN!(data, data_type, solution, stats)

                finish_HC_time = time_ns()
                total_time_hc = (finish_HC_time - start_time) * 1e-9
                stats.HC_time = total_time_hc

            elseif params.approach_HC == "IMB"

                start_time = time_ns()

                # Setup std formulation solution data structure
                solution = OutputStatistics.init_std_form_solution(data_type)

                #solver for VMP

                Constructive.solve_IMB!(data, data_type, solution, stats)

                finish_HC_time = time_ns()
                total_time_hc = (finish_HC_time - start_time) * 1e-9
                stats.HC_time = total_time_hc
            elseif params.approach_HC == "random"
                start_time = time_ns()

                # Setup std formulation solution data structure
                solution = OutputStatistics.init_std_form_solution(data_type)

                #solver for VMP

                Constructive.solve_random!(data, data_type, solution, stats)

                finish_HC_time = time_ns()
                total_time_hc = (finish_HC_time - start_time) * 1e-9
                stats.HC_time = total_time_hc
            elseif params.approach_HC == "C_GRASP"
                start_time = time_ns()

                # Setup std formulation solution data structure
                solution = OutputStatistics.init_std_form_solution(data_type)

                #solver for VMP

                Constructive.solve_C_GRASP!(data, data_type, solution, stats)

                finish_HC_time = time_ns()
                total_time_hc = (finish_HC_time - start_time) * 1e-9
                stats.HC_time = total_time_hc
            end


            
            if params.approach == "HC_LS"

                if params.approach_LS == "bestImprovement"
                  


                    if params.approach_improvemente == "Swap"
                        start_time_ls = time_ns()
                        #solver for VMP
                        BestImprovement.bestImprovementSwap!(data, data_type, solution, stats)
                        finish_LS_time = time_ns()
                        total_time_ls = (finish_LS_time - start_time_ls) * 1e-9
                        stats.LS_time = total_time_ls
                    elseif params.approach_improvemente == "Reinsertion"
                        start_time_ls = time_ns()
                        #solver for VMP
                        BestImprovement.bestImprovementReinsertion!(data, data_type, solution, stats)
                        finish_LS_time = time_ns()
                        total_time_ls = (finish_LS_time - start_time_ls) * 1e-9
                        stats.LS_time = total_time_ls
                    elseif params.approach_improvemente == "Or-opt-2"
                        start_time_ls = time_ns()
                        #solver for VMP
                        BestImprovement.bestImprovementOrOpt2!(data, data_type, solution, stats)
                        finish_LS_time = time_ns()
                        total_time_ls = (finish_LS_time - start_time_ls) * 1e-9
                        stats.LS_time = total_time_ls
                    elseif params.approach_improvemente == "Or-opt-3"
                        start_time_ls = time_ns()
                        #solver for VMP
                        BestImprovement.bestImprovementOrOpt3!(data, data_type, solution, stats)
                        finish_LS_time = time_ns()
                        total_time_ls = (finish_LS_time - start_time_ls) * 1e-9
                        stats.LS_time = total_time_ls
                    elseif params.approach_improvemente == "2-opt"
                        start_time_ls = time_ns()
                        #solver for VMP
                        BestImprovement.bestImprovement2Opt!(data, data_type, solution, stats)
                        finish_LS_time = time_ns()
                        total_time_ls = (finish_LS_time - start_time_ls) * 1e-9
                        stats.LS_time = total_time_ls
                        
                    end
                elseif params.approach_LS ==  "RVND"
                    start_time_ls = time_ns()
                    #solver for VMP
                    maxIter = 10
                    Iter = 0
                    while Iter < maxIter
                        
                        
                        RD = rand(1:5)
                        RD = 5  
                        if RD == 1
                            BestImprovement.bestImprovementSwap!(data, data_type, solution, stats)
                        elseif RD == 2
                            BestImprovement.bestImprovementReinsertion!(data, data_type, solution, stats)
                        elseif RD == 3
                            BestImprovement.bestImprovementOrOpt2!(data, data_type, solution, stats)
                        elseif RD == 4
                           BestImprovement.bestImprovementOrOpt3!(data, data_type, solution, stats)
                        elseif RD == 5
                            BestImprovement.bestImprovement2Opt!(data, data_type, solution, stats)
                        end

                        if stats.improvement == 0
                            Iter += 1
                        else
                            Iter = 0
                        end
                        #println("stats.UB_LS aqui !!!!", stats.UB_LS)

                        if stats.UB_LS !=  -1e8
                            if stats.UB_LS < 0 #12
                                println("parar")
                                break
                            end 
                        end
                        

                    end
                    
                    
                    
                    
                    stats.improvement = ((stats.UB_HC - stats.UB_LS)/stats.UB_HC)*100

                    finish_LS_time = time_ns()
                    total_time_ls = (finish_LS_time - start_time_ls) * 1e-9
                    stats.LS_time = total_time_ls

                end
            end
 
        elseif params.approach == "ILS"   

            start_time_ls = start_time = time_ns()

            # Setup std formulation solution data structure
            solution = OutputStatistics.init_std_form_solution(data_type)

            #solver for VMP

            MetaHeuristics.ILS!(data, data_type, solution, stats)

            stats.improvement = ((stats.UB_HC - stats.UB_LS)/stats.UB_HC)*100

            finish_LS_time = time_ns()
            total_time_ls = (finish_LS_time - start_time_ls) * 1e-9
            stats.LS_time = total_time_ls
        elseif params.approach == "MultiStart"   

            start_time_ls = start_time = time_ns()

            # Setup std formulation solution data structure
            solution = OutputStatistics.init_std_form_solution(data_type)

            #solver for VMP

            MetaHeuristics.MultiStart!(data, data_type, solution, stats)

            stats.improvement = ((stats.UB_HC - stats.UB_LS)/stats.UB_HC)*100

            finish_LS_time = time_ns()
            total_time_ls = (finish_LS_time - start_time_ls) * 1e-9
            stats.LS_time = total_time_ls  
        elseif params.approach == "SA"   

            start_time_ls = start_time = time_ns()

            # Setup std formulation solution data structure
            solution = OutputStatistics.init_std_form_solution(data_type)

            #solver for VMP

            MetaHeuristics.SA!(data, data_type, solution, stats)

            stats.improvement = ((stats.UB_HC - stats.UB_LS)/stats.UB_HC)*100

            finish_LS_time = time_ns()
            total_time_ls = (finish_LS_time - start_time_ls) * 1e-9
            stats.LS_time = total_time_ls  
        elseif params.approach == "VNS"   

            start_time_ls = start_time = time_ns()

            # Setup std formulation solution data structure
            solution = OutputStatistics.init_std_form_solution(data_type)

            #solver for VMP

            MetaHeuristics.VNS!(data, data_type, solution, stats)

            stats.improvement = ((stats.UB_HC - stats.UB_LS)/stats.UB_HC)*100

            finish_LS_time = time_ns()
            total_time_ls = (finish_LS_time - start_time_ls) * 1e-9
            stats.LS_time = total_time_ls    
        elseif params.approach == "GRASP"   

            start_time_ls = start_time = time_ns()

            # Setup std formulation solution data structure
            solution = OutputStatistics.init_std_form_solution(data_type)

            #solver for VMP

            MetaHeuristics.GRASP!(data, data_type, solution, stats)

            stats.improvement = ((stats.UB_HC - stats.UB_LS)/stats.UB_HC)*100

            finish_LS_time = time_ns()
            total_time_ls = (finish_LS_time - start_time_ls) * 1e-9
            stats.LS_time = total_time_ls 
        elseif params.approach == "TABU"   

            start_time_ls = start_time = time_ns()

            # Setup std formulation solution data structure
            solution = OutputStatistics.init_std_form_solution(data_type)

            #solver for VMP

            MetaHeuristics.TABU!(data, data_type, solution, stats)

            stats.improvement = ((stats.UB_HC - stats.UB_LS)/stats.UB_HC)*100

            finish_LS_time = time_ns()
            total_time_ls = (finish_LS_time - start_time_ls) * 1e-9
            stats.LS_time = total_time_ls
        end   
 

        finish_time = time_ns()
        total_time = (finish_time - start_time) * 1e-9
        stats.total_time = total_time
        # Print statistics to output file
        OutputStatistics.print_stats(data_type, params, stats, output_file)
        
        
    

    end


end

main(ARGS)