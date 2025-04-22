module FixAndOptimize

using JuMP
using Gurobi
using Data
using ReformulationFL
using ReformulationSR

using Formulation
using OutputStatistics
using Parameters
using CPUTime


export fix_and_optimize_per_period!, fix_and_optimize_per_period2!, fix_and_optimize_per_period_per_item!, fix_and_optimize_per_period_0!, fix_and_optimize_per_period_1!, fix_and_optimize_per_period_Relaxation_Search!, fix_and_optimize_per_period_Relaxation_Search2

function fix_and_optimize_per_period!(data::InstanceData, model::Model,
    sol::Union{StdFormModelSolution, FlFormModelSolution, SrFormModelSolution}, params::ExperimentParameters, stats::StatisticsData)

    println("===================== Running FO... =====================")

    I = 1:data.NI
    T = 1:data.NT
    

    elapsed_time = 0.0

    # Get fix-and-optimize start time
    time_start = time_ns()

    stats.fo_iterations = 0
    stats.fo_improvement = 0.0
    stats.fo_time = 0.0

    stop_by_time_limit = false

    best_UB = sol.primal_bound
    println("Current primal bound = ", best_UB)
    curr_UB = best_UB
    improved = true
    improvement = 0.0

    # Inititalize window
    window_size = Int64(max(2, round(data.NT / params.fo_window_size_factor)))
    step_size = Int64(round(window_size / params.fo_step_size_factor))
    max_window_size = Int64(round(data.NT / params.fo_max_window_size_factor))#  params.fo_max_window_size_factor * data.NT




    if params.fo_theta == 0
        params.fo_max_rounds_without_improvement = 1
    end

    # sol_is_feasible = Formulation.check_solution(data, sol)

    # if sol_is_feasible == true
    #     println("Solution is FEASIBLE")
    #     for k in M, t in T
    #         set_upper_bound(model[:cap_viol_var][t], 0.0)
    #     end
    # else
    #     println("Solution is INFEASIBLE")
    # end

    not_improved_rounds = 0
    stats.fo_rounds = 0
    params.fo_max_time = params.total_time_limit - stats.rf_time
    available_time = params.fo_max_time  
    println("available time: ", available_time, " RF solution = ", best_UB)

    Parameters.set_fo_parameters(model, params)

    while elapsed_time < params.fo_max_time && window_size <= max_window_size && not_improved_rounds < params.fo_max_rounds_without_improvement
        # println("rodar a primeira parte")
        stats.fo_rounds += 1

        # Inititalize window
        window_begin = 1

        # Compute end of window
        window_end = window_begin + window_size - 1

        improved = false
        improvement = 0.0

        #=  for t in window_begin:window_end
            for i in I
                println("Upper bound of w[$i, $t] = ", upper_bound(model[:w][i, t]))
            end
        end =#

        while window_end <= data.NT && window_begin <= data.NT && elapsed_time < params.fo_max_time  && sol.status != INFEASIBLE

            stats.fo_iterations += 1
            print( "\nrnd = ", stats.fo_rounds, " iter: ",  stats.fo_iterations,  " [$window_begin..$window_end]", " best_UB = ", round(best_UB, digits = 5))

            # Fix variables before the window
            if window_begin > 1
                for t in 1:window_begin - 1
                    for i in I
                        
                        if sol.y[i, t] >= 0.9
                            set_lower_bound(model[:y][i, t], 1.0)
                            set_upper_bound(model[:y][i, t], 1.0)
                        end
                        if sol.y[i, t] <= 0.1
                            set_lower_bound(model[:y][i, t], 0.0)
                            set_upper_bound(model[:y][i, t], 0.0)
                        end
                        if sol.w[i, t] >= 0.9
                            set_lower_bound(model[:w][i, t], 1.0)
                            set_upper_bound(model[:w][i, t], 1.0)
                        end
                        if sol.w[i, t] <= 0.1
                            set_lower_bound(model[:w][i, t], 0.0)
                            set_upper_bound(model[:w][i, t], 0.0)
                        end
                        
                    end
                end
            end

            # Define variables as binary inside the window
            for t in window_begin:window_end
                for i in I
                    
                    if !is_binary(model[:y][i, t])
                        set_binary(model[:y][i, t])
                    end
                    set_lower_bound(model[:y][i, t], 0.0)
                    set_upper_bound(model[:y][i, t], 1.0)

                    if !is_binary(model[:w][i, t])
                        set_binary(model[:w][i, t])
                    end
                    set_lower_bound(model[:w][i, t], 0.0)
                    set_upper_bound(model[:w][i, t], 1.0)

                    # set_start_value(model[:y][i,k,t], sol.y[i,k,t])
                    
                end
            end

            # Fix variables after the window
            if window_begin == 1
                for t in window_end + 1:data.NT
                    for i in I

                        if sol.y[i, t] >= 0.9
                            set_lower_bound(model[:y][i, t], 1.0)
                            set_upper_bound(model[:y][i, t], 1.0)
                        end
                        if sol.y[i, t] <= 0.1
                            set_lower_bound(model[:y][i, t], 0.0)
                            set_upper_bound(model[:y][i, t], 0.0)
                        end
                        if sol.w[i, t] >= 0.9
                            set_lower_bound(model[:w][i, t], 1.0)
                            set_upper_bound(model[:w][i, t], 1.0)
                        end
                        if sol.w[i, t] <= 0.1
                            set_lower_bound(model[:w][i, t], 0.0)
                            set_upper_bound(model[:w][i, t], 0.0)
                        end

                    end
                end
            end

            # Solve restricted model
            set_time_limit_sec(model, min(available_time, params.fo_restricted_model_max_time))
            #write_to_file(model, "model.lp")

            if params.MIP_model == "std_form"
                println("Solving std_form model")
                Formulation.solve_stdform_model!(model, data, sol, stats)
                
            elseif params.MIP_model == "fac_loc_form"
                println("Solving fl_form model")
                ReformulationFL.solve_flform_model!(model, data, sol, stats)
            elseif params.MIP_model == "short_path_form"
                println("Solving rs_form model")
                ReformulationSR.solve_rsform_model!(model, data, sol, stats)
            end
             
            sol.status = termination_status(model)
            println(sol.status)
            curr_UB = sol.primal_bound
            improvement = 0.0

            if best_UB - curr_UB > 0.000001
                improved = true
                not_improved_rounds = 0
                improvement = 100 * ((best_UB - curr_UB) / best_UB)
                best_UB = curr_UB
            end

            print(" curr_obj = ", round(curr_UB, digits = 5),  " improvement = ", round(improvement, digits = 4))

            curr_time = time_ns()
            elapsed_time = ((curr_time - time_start) * 1e-9)
            print(" elapsed_time = ", round(elapsed_time, digits = 2))

            available_time = params.fo_max_time - elapsed_time
            print(" available time: ", round(available_time, digits = 4))

            if elapsed_time > params.fo_max_time
                stop_by_time_limit = true
                println("\nFO time limit exceeded")
                break
            end
            
            if window_end == data.NT
                break
            end

            # Update FO window
            window_begin += step_size
            window_end = min(window_begin + window_size - 1, data.NT)

        end

        println("\nSolution was improved in round $(stats.fo_rounds)? R.: ", improved)
        curr_time = time_ns()
        elapsed_time = ((curr_time - time_start) * 1e-9)

        available_time = params.fo_max_time - elapsed_time

        if elapsed_time > params.fo_max_time
            stop_by_time_limit = true
            println("\nFO time limit exceeded")
            break
        end

        # Compare not improved rounds for theta = 0
        if !improved && params.fo_theta == 0
            not_improved_rounds += 1
        end

        # Dynamically changes window_size and step_size        
        if !improved && params.fo_theta > 0
            not_improved_rounds += 1
            window_size += params.fo_theta
            if window_size - step_size > 2
                step_size += 1
            end
        else # Reset window_size and step_size
            window_size = Int64(max(2, floor(data.NT / params.fo_window_size_factor)))
            step_size = Int64(ceil(window_size / params.fo_step_size_factor))
        end
    end

    # Check if solution is feasible (integer)
    integerSol = true

    for i in I, t in T
        if value(sol.y[i, t]) > 0.01 && value(sol.y[i, t]) < 0.99
            integerSol = false
            println("FFAAIILL!!")
            break
        end
    end

    if integerSol == true
        if stop_by_time_limit == false
            stats.sol_status = "Success"
            println("\nRF found an integer solution with value: ", sol.primal_bound)
        else
            stats.sol_status = "Success_TL"
            println("\nRF stopped by time limit but found an integer solution with value: ", sol.primal_bound)
        end
    else
        if stop_by_time_limit == false
            stats.sol_status = "Fail: !Int"
            println("\nRF could not find an integer solution")
        else
            stats.sol_status = "Fail: !Int_TL"
            println("\nRF stopped by time limit and could not find an integer solution")
        end
    end

    time_end = time_ns()

    stats.fo_time = ((time_end - time_start) * 1e-9)
    stats.fo_UB = best_UB
    stats.fo_improvement = 100 * ((stats.rf_UB - stats.fo_UB) / stats.rf_UB)

    println("\nBest primal bound = ", best_UB)

    println("\nfo_time = ", stats.fo_time)

    println("===================== Finished FO! ======================")

end
function fix_and_optimize_per_period2!(data::InstanceData, model::Model,
    sol::Union{StdFormModelSolution, FlFormModelSolution, SrFormModelSolution}, params::ExperimentParameters, stats::StatisticsData)

    println("===================== Running FO... =====================")

    I = 1:data.NI
    T = 1:data.NT
    

    elapsed_time = 0.0

    # Get fix-and-optimize start time
    time_start = time_ns()

    stats.fo_iterations = 0
    stats.fo_improvement = 0.0
    stats.fo_time = 0.0

    stop_by_time_limit = false

    best_UB = sol.primal_bound
    println("Current primal bound = ", best_UB)
    curr_UB = best_UB
    improved = true
    improvement = 0.0

    # Inititalize window
    window_size = Int64(max(2, round(data.NT / params.fo_window_size_factor)))
    step_size = Int64(round(window_size / params.fo_step_size_factor))
    max_window_size = Int64(round(data.NT / params.fo_max_window_size_factor))#  params.fo_max_window_size_factor * data.NT
    q_mais = 1
    window_size = window_size + q_mais



    if params.fo_theta == 0
        params.fo_max_rounds_without_improvement = 1
    end

    # sol_is_feasible = Formulation.check_solution(data, sol)

    # if sol_is_feasible == true
    #     println("Solution is FEASIBLE")
    #     for k in M, t in T
    #         set_upper_bound(model[:cap_viol_var][t], 0.0)
    #     end
    # else
    #     println("Solution is INFEASIBLE")
    # end

    not_improved_rounds = 0
    stats.fo_rounds = 0
    params.fo_max_time = params.total_time_limit - stats.rf_time
    available_time = params.fo_max_time  
    println("available time: ", available_time, " RF solution = ", best_UB)

    Parameters.set_fo_parameters(model, params)

    while elapsed_time < params.fo_max_time && window_size <= max_window_size && not_improved_rounds < params.fo_max_rounds_without_improvement
        # println("rodar a primeira parte")
        stats.fo_rounds += 1

        # Inititalize window
        window_begin = 1
        window_begin2 = max(window_begin - q_mais, 1)
        # Compute end of window
        window_end = window_begin + window_size - 1
        window_end2 = min(window_end + q_mais, data.NT)
        improved = false
        improvement = 0.0

        #=  for t in window_begin:window_end
            for i in I
                println("Upper bound of w[$i, $t] = ", upper_bound(model[:w][i, t]))
            end
        end =#

        while window_end <= data.NT && window_begin <= data.NT && elapsed_time < params.fo_max_time  && sol.status != INFEASIBLE

            stats.fo_iterations += 1
            print( "\nrnd = ", stats.fo_rounds, " iter: ",  stats.fo_iterations,  " [$window_begin..$window_end]", " best_UB = ", round(best_UB, digits = 5))

            # Fix variables before the window
            if window_begin > 1
                for t in 1:window_begin - 1
                    for i in I
                        
                        if sol.y[i, t] >= 0.9
                            set_lower_bound(model[:y][i, t], 1.0)
                            set_upper_bound(model[:y][i, t], 1.0)
                        end
                        if sol.y[i, t] <= 0.1
                            set_lower_bound(model[:y][i, t], 0.0)
                            set_upper_bound(model[:y][i, t], 0.0)
                        end

                        
                    end
                end
            end
            if window_begin2 > 1 
                for t in 1:window_begin - 1
                    for i in I
                        
                        
                        if sol.w[i, t] >= 0.9
                            set_lower_bound(model[:w][i, t], 1.0)
                            set_upper_bound(model[:w][i, t], 1.0)
                        end
                        if sol.w[i, t] <= 0.1
                            set_lower_bound(model[:w][i, t], 0.0)
                            set_upper_bound(model[:w][i, t], 0.0)
                        end
                        
                    end
                end
            end

            # Define variables as binary inside the window
            for t in window_begin:window_end
                for i in I
                    
                    if !is_binary(model[:y][i, t])
                        set_binary(model[:y][i, t])
                    end
                    set_lower_bound(model[:y][i, t], 0.0)
                    set_upper_bound(model[:y][i, t], 1.0)

      
                   
                    
                end
            end
            for t in window_begin2:window_end2
                for i in I
                    


                    if !is_binary(model[:w][i, t])
                        set_binary(model[:w][i, t])
                    end
                    set_lower_bound(model[:w][i, t], 0.0)
                    set_upper_bound(model[:w][i, t], 1.0)

                   
                    
                end
            end

            # Fix variables after the window
            if window_begin == 1
                for t in window_end + 1:data.NT
                    for i in I

                        if sol.y[i, t] >= 0.9
                            set_lower_bound(model[:y][i, t], 1.0)
                            set_upper_bound(model[:y][i, t], 1.0)
                        end
                        if sol.y[i, t] <= 0.1
                            set_lower_bound(model[:y][i, t], 0.0)
                            set_upper_bound(model[:y][i, t], 0.0)
                        end
                        if sol.w[i, t] >= 0.9
                            set_lower_bound(model[:w][i, t], 1.0)
                            set_upper_bound(model[:w][i, t], 1.0)
                        end
                        if sol.w[i, t] <= 0.1
                            set_lower_bound(model[:w][i, t], 0.0)
                            set_upper_bound(model[:w][i, t], 0.0)
                        end

                    end
                end
            end

            # Solve restricted model
            set_time_limit_sec(model, min(available_time, params.fo_restricted_model_max_time))
            #write_to_file(model, "model.lp")

            if params.MIP_model == "std_form"
                println("Solving std_form model")
                Formulation.solve_stdform_model!(model, data, sol, stats)
                
            elseif params.MIP_model == "fac_loc_form"
                println("Solving fl_form model")
                ReformulationFL.solve_flform_model!(model, data, sol, stats)
            elseif params.MIP_model == "short_path_form"
                println("Solving rs_form model")
                ReformulationSR.solve_rsform_model!(model, data, sol, stats)
            end
             
            sol.status = termination_status(model)
            println(sol.status)
            curr_UB = sol.primal_bound
            improvement = 0.0

            if best_UB - curr_UB > 0.000001
                improved = true
                not_improved_rounds = 0
                improvement = 100 * ((best_UB - curr_UB) / best_UB)
                best_UB = curr_UB
            end

            print(" curr_obj = ", round(curr_UB, digits = 5),  " improvement = ", round(improvement, digits = 4))

            curr_time = time_ns()
            elapsed_time = ((curr_time - time_start) * 1e-9)
            print(" elapsed_time = ", round(elapsed_time, digits = 2))

            available_time = params.fo_max_time - elapsed_time
            print(" available time: ", round(available_time, digits = 4))

            if elapsed_time > params.fo_max_time
                stop_by_time_limit = true
                println("\nFO time limit exceeded")
                break
            end
            
            if window_end == data.NT
                break
            end

            # Update FO window
            window_begin += step_size
            window_end = min(window_begin + window_size - 1, data.NT)
            window_begin2 = min(window_begin - q_mais, 1)
            window_end2 = min(window_end2 + q_mais, data.NT)

            

        end

        println("\nSolution was improved in round $(stats.fo_rounds)? R.: ", improved)
        curr_time = time_ns()
        elapsed_time = ((curr_time - time_start) * 1e-9)

        available_time = params.fo_max_time - elapsed_time

        if elapsed_time > params.fo_max_time
            stop_by_time_limit = true
            println("\nFO time limit exceeded")
            break
        end

        # Compare not improved rounds for theta = 0
        if !improved && params.fo_theta == 0
            not_improved_rounds += 1
        end

        # Dynamically changes window_size and step_size        
        if !improved && params.fo_theta > 0
            not_improved_rounds += 1
            window_size += params.fo_theta
            if window_size - step_size > 2
                step_size += 1
            end
        else # Reset window_size and step_size
            window_size = Int64(max(2, floor(data.NT / params.fo_window_size_factor)))
            step_size = Int64(ceil(window_size / params.fo_step_size_factor))
        end
    end

    # Check if solution is feasible (integer)
    integerSol = true

    for i in I, t in T
        if value(sol.y[i, t]) > 0.01 && value(sol.y[i, t]) < 0.99
            integerSol = false
            println("FFAAIILL!!")
            break
        end
    end

    if integerSol == true
        if stop_by_time_limit == false
            stats.sol_status = "Success"
            println("\nRF found an integer solution with value: ", sol.primal_bound)
        else
            stats.sol_status = "Success_TL"
            println("\nRF stopped by time limit but found an integer solution with value: ", sol.primal_bound)
        end
    else
        if stop_by_time_limit == false
            stats.sol_status = "Fail: !Int"
            println("\nRF could not find an integer solution")
        else
            stats.sol_status = "Fail: !Int_TL"
            println("\nRF stopped by time limit and could not find an integer solution")
        end
    end

    time_end = time_ns()

    stats.fo_time = ((time_end - time_start) * 1e-9)
    stats.fo_UB = best_UB
    stats.fo_improvement = 100 * ((stats.rf_UB - stats.fo_UB) / stats.rf_UB)

    println("\nBest primal bound = ", best_UB)

    println("\nfo_time = ", stats.fo_time)

    println("===================== Finished FO! ======================")

end
function fix_and_optimize_per_period_per_item!(data::InstanceData, model::Model,
    sol::Union{StdFormModelSolution, FlFormModelSolution, SrFormModelSolution}, params::ExperimentParameters, stats::StatisticsData)

    println("===================== Running FO... =====================")

    I = 1:data.NI
    T = 1:data.NT
    

    elapsed_time = 0.0

    # Get fix-and-optimize start time
    time_start = time_ns()

    stats.fo_iterations = 0
    stats.fo_improvement = 0.0
    stats.fo_time = 0.0

    stop_by_time_limit = false

    best_UB = sol.primal_bound
    println("Current primal bound = ", best_UB)
    curr_UB = best_UB
    improved = true
    improvement = 0.0

    # Inititalize window
    window_size = Int64(max(2, round(data.NT / params.fo_window_size_factor)))
    step_size = Int64(round(window_size / params.fo_step_size_factor))
    max_window_size = Int64(round(data.NT / params.fo_max_window_size_factor))#  params.fo_max_window_size_factor * data.NT




    if params.fo_theta == 0
        params.fo_max_rounds_without_improvement = 1
    end

    # sol_is_feasible = Formulation.check_solution(data, sol)

    # if sol_is_feasible == true
    #     println("Solution is FEASIBLE")
    #     for k in M, t in T
    #         set_upper_bound(model[:cap_viol_var][t], 0.0)
    #     end
    # else
    #     println("Solution is INFEASIBLE")
    # end

    not_improved_rounds = 0

    params.fo_max_time = params.total_time_limit - stats.rf_time
    available_time = params.fo_max_time  
    println("available time: ", available_time, " RF solution = ", best_UB)

    Parameters.set_fo_parameters(model, params)

    while elapsed_time < params.fo_max_time && window_size <= max_window_size && not_improved_rounds < params.fo_max_rounds_without_improvement
        println("rodar a primeira parte")
        stats.fo_rounds += 1

        # Inititalize window
        window_begin = 1

        # Compute end of window
        window_end = window_begin + window_size - 1

        improved = false
        improvement = 0.0



        while window_end <= data.NT && window_begin <= data.NT && elapsed_time < params.fo_max_time  && sol.status != INFEASIBLE

            stats.fo_iterations += 1
            print( "\nrnd = ", stats.fo_rounds, " iter: ",  stats.fo_iterations,  " [$window_begin..$window_end]", " best_UB = ", round(best_UB, digits = 5))

            # Fix variables before the window
            if window_begin > 1
                for t in 1:window_begin - 1
                    for i in I
                        
                        if sol.y[i, t] >= 0.9
                            set_lower_bound(model[:y][i, t], 1.0)
                        end
                        if sol.y[i, t] <= 0.1
                            set_upper_bound(model[:y][i, t], 0.0)
                        end
                        if sol.w[i, t] >= 0.9
                            set_lower_bound(model[:w][i, t], 1.0)
                        end
                        if sol.w[i, t] <= 0.1
                            set_upper_bound(model[:w][i, t], 0.0)
                        end
                        
                    end
                end
            end

            # Define variables as binary inside the window
            for t in window_begin:window_end
                for i in I
                    
                    if !is_binary(model[:y][i, t])
                        set_binary(model[:y][i, t])
                    end
                    set_lower_bound(model[:y][i, t], 0.0)
                    set_upper_bound(model[:y][i, t], 1.0)

                    if !is_binary(model[:w][i, t])
                        set_binary(model[:w][i, t])
                    end
                    set_lower_bound(model[:w][i, t], 0.0)
                    set_upper_bound(model[:w][i, t], 1.0)

                    # set_start_value(model[:y][i,k,t], sol.y[i,k,t])
                    
                end
            end

            # Fix variables after the window
            if window_begin == 1
                for t in window_end + 1:data.NT
                    for i in I

                        if sol.y[i, t] >= 0.9
                            set_lower_bound(model[:y][i, t], 1.0)
                        end
                        if sol.y[i, t] <= 0.1
                            set_upper_bound(model[:y][i, t], 0.0)
                        end
                        if sol.w[i, t] >= 0.9
                            set_lower_bound(model[:w][i, t], 1.0)
                        end
                        if sol.w[i, t] <= 0.1
                            set_upper_bound(model[:w][i, t], 0.0)
                        end

                    end
                end
            end

            # Solve restricted model
            set_time_limit_sec(model, min(available_time, params.fo_restricted_model_max_time))
            #write_to_file(model, "model.lp")

            if params.MIP_model == "std_form"
                println("Solving std_form model")
                Formulation.solve_stdform_model!(model, data, sol, stats)
            elseif params.MIP_model == "fac_loc_form"
                println("Solving fl_form model")
                ReformulationFL.solve_flform_model!(model, data, sol, stats)
            elseif params.MIP_model == "short_path_form"
                println("Solving rs_form model")
                ReformulationSR.solve_rsform_model!(model, data, sol, stats)
            end
             
            sol.status = termination_status(model)
            println(sol.status)
            curr_UB = sol.primal_bound
            improvement = 0.0

            if best_UB - curr_UB > 0.000001
                improved = true
                not_improved_rounds = 0
                improvement = 100 * ((best_UB - curr_UB) / best_UB)
                best_UB = curr_UB
            end

            print(" curr_obj = ", round(curr_UB, digits = 5),  " improvement = ", round(improvement, digits = 4))

            curr_time = time_ns()
            elapsed_time = ((curr_time - time_start) * 1e-9)
            print(" elapsed_time = ", round(elapsed_time, digits = 2))

            available_time = params.fo_max_time - elapsed_time
            print(" available time: ", round(available_time, digits = 4))

            if elapsed_time > params.fo_max_time
                stop_by_time_limit = true
                println("\nFO time limit exceeded")
                break
            end
            
            if window_end == data.NT
                break
            end

            # Update FO window
            window_begin += step_size
            window_end = min(window_begin + window_size - 1, data.NT)

        end

        println("\nSolution was improved in round $(stats.fo_rounds)? R.: ", improved)
        curr_time = time_ns()
        elapsed_time = ((curr_time - time_start) * 1e-9)

        available_time = params.fo_max_time - elapsed_time

        if elapsed_time > params.fo_max_time
            stop_by_time_limit = true
            println("\nFO time limit exceeded")
            break
        end

        # Compare not improved rounds for theta = 0
        if !improved && params.fo_theta == 0
            not_improved_rounds += 1
        end

        # Dynamically changes window_size and step_size        
        if !improved && params.fo_theta > 0
            not_improved_rounds += 1
            window_size += params.fo_theta
            if window_size - step_size > 2
                step_size += 1
            end
        else # Reset window_size and step_size
            window_size = Int64(max(2, floor(data.NT / params.fo_window_size_factor)))
            step_size = Int64(ceil(window_size / params.fo_step_size_factor))
        end
    end

    # Check if solution is feasible (integer)
    integerSol = true

    for i in I, t in T
        if value(sol.y[i, t]) > 0.01 && value(sol.y[i, t]) < 0.99
            integerSol = false
            println("FFAAIILL!!")
            break
        end
    end

    if integerSol == true
        if stop_by_time_limit == false
            stats.sol_status = "Success"
            println("\nRF found an integer solution with value: ", sol.primal_bound)
        else
            stats.sol_status = "Success_TL"
            println("\nRF stopped by time limit but found an integer solution with value: ", sol.primal_bound)
        end
    else
        if stop_by_time_limit == false
            stats.sol_status = "Fail: !Int"
            println("\nRF could not find an integer solution")
        else
            stats.sol_status = "Fail: !Int_TL"
            println("\nRF stopped by time limit and could not find an integer solution")
        end
    end

    time_end = time_ns()

    stats.fo_time = ((time_end - time_start) * 1e-9)
    stats.fo_UB = best_UB
    stats.fo_improvement = 100 * ((stats.rf_UB - stats.fo_UB) / stats.rf_UB)

    println("\nBest primal bound = ", best_UB)

    println("\nfo_time = ", stats.fo_time)

    println("===================== Finished FO! ======================")

end
function fix_and_optimize_per_period_0!(data::InstanceData, model::Model,
    sol::Union{StdFormModelSolution, FlFormModelSolution, SrFormModelSolution}, params::ExperimentParameters, stats::StatisticsData)

    println("===================== Running FO... =====================")

    I = 1:data.NI
    T = 1:data.NT
    

    elapsed_time = 0.0

    # Get fix-and-optimize start time
    time_start = time_ns()

    stats.fo_iterations = 0
    stats.fo_improvement = 0.0
    stats.fo_time = 0.0

    stop_by_time_limit = false

    best_UB = sol.primal_bound
    println("Current primal bound = ", best_UB)
    curr_UB = best_UB
    improved = true
    improvement = 0.0

    # Inititalize window
    window_size = Int64(max(2, round(data.NT / params.fo_window_size_factor2)))
    step_size = Int64(round(window_size / params.fo_step_size_factor2))
    max_window_size = Int64(round(data.NT / params.fo_max_window_size_factor2))#  params.fo_max_window_size_factor * data.NT




    if params.fo_theta2 == 0
        params.fo_max_rounds_without_improvement2 = 1
    end

    # sol_is_feasible = Formulation.check_solution(data, sol)

    # if sol_is_feasible == true
    #     println("Solution is FEASIBLE")
    #     for k in M, t in T
    #         set_upper_bound(model[:cap_viol_var][t], 0.0)
    #     end
    # else
    #     println("Solution is INFEASIBLE")
    # end

    not_improved_rounds = 0

    params.fo_max_time = params.total_time_limit - stats.rf_time
    available_time = params.fo_max_time  
    println("available time: ", available_time, " RF solution = ", best_UB)

    Parameters.set_fo_parameters(model, params)




    while elapsed_time < params.fo_max_time && window_size <= max_window_size && not_improved_rounds < params.fo_max_rounds_without_improvement2
        println("rodar a primeira parte")
        stats.fo_rounds += 1

        # Inititalize window
        window_begin = 1

        # Compute end of window
        window_end = window_begin + window_size - 1

        improved = false
        improvement = 0.0

        for i in I
            for t in T
                
                set_upper_bound(model[:w][i, t], 0.0)
                
            end
        end


        while window_end <= data.NT && window_begin <= data.NT && elapsed_time < params.fo_max_time  && sol.status != INFEASIBLE

            stats.fo_iterations += 1
            print( "\nrnd = ", stats.fo_rounds, " iter: ",  stats.fo_iterations,  " [$window_begin..$window_end]", " best_UB = ", round(best_UB, digits = 5))

            # Fix variables before the window
            if window_begin > 1
                for t in 1:window_begin - 1
                    for i in I
                        
                        if sol.y[i, t] >= 0.9
                            set_lower_bound(model[:y][i, t], 1.0)
                            set_upper_bound(model[:y][i, t], 1.0)
                        end
                        if sol.y[i, t] <= 0.1
                            set_lower_bound(model[:y][i, t], 0.0)
                            set_upper_bound(model[:y][i, t], 0.0)
                        end
                        
                    end
                end
            end

            # Define variables as binary inside the window
            for t in window_begin:window_end
                for i in I
                    
                    if !is_binary(model[:y][i, t])
                        set_binary(model[:y][i, t])
                    end
                    set_lower_bound(model[:y][i, t], 0.0)
                    set_upper_bound(model[:y][i, t], 1.0)

                    # set_start_value(model[:y][i,k,t], sol.y[i,k,t])
                    
                end
            end

            # Fix variables after the window
            if window_begin == 1
                for t in window_end + 1:data.NT
                    for i in I

                        if sol.y[i, t] >= 0.9
                            set_lower_bound(model[:y][i, t], 1.0)
                            set_upper_bound(model[:y][i, t], 1.0)
                        end
                        if sol.y[i, t] <= 0.1
                            set_lower_bound(model[:y][i, t], 0.0)
                            set_upper_bound(model[:y][i, t], 0.0)
                        end

                    end
                end
            end

            # Solve restricted model
            set_time_limit_sec(model, min(available_time, params.fo_restricted_model_max_time))
            #write_to_file(model, "model.lp")

            if params.MIP_model == "std_form"
                println("Solving std_form model")
                Formulation.solve_stdform_model!(model, data, sol, stats)
            elseif params.MIP_model == "fac_loc_form"
                println("Solving fl_form model")
                ReformulationFL.solve_flform_model!(model, data, sol, stats)
            elseif params.MIP_model == "short_path_form"
                println("Solving rs_form model")
                ReformulationSR.solve_rsform_model!(model, data, sol, stats)
            end
             
            sol.status = termination_status(model)
            println(sol.status)
            curr_UB = sol.primal_bound
            improvement = 0.0

            #= if best_UB - curr_UB > 0.000001
                improved = true
                not_improved_rounds = 0
                improvement = 100 * ((best_UB - curr_UB) / best_UB)
                best_UB = curr_UB
            end =#

            print(" curr_obj = ", round(curr_UB, digits = 5),  " improvement = ", round(improvement, digits = 4))

            curr_time = time_ns()
            elapsed_time = ((curr_time - time_start) * 1e-9)
            print(" elapsed_time = ", round(elapsed_time, digits = 2))

            available_time = params.fo_max_time - elapsed_time
            print(" available time: ", round(available_time, digits = 4))

            if elapsed_time > params.fo_max_time
                stop_by_time_limit = true
                println("\nFO time limit exceeded")
                break
            end
            
            if window_end == data.NT
                break
            end

            # Update FO window
            window_begin += step_size
            window_end = min(window_begin + window_size - 1, data.NT)

        end

        println("\nSolution was improved in round $(stats.fo_rounds)? R.: ", improved)
        curr_time = time_ns()
        elapsed_time = ((curr_time - time_start) * 1e-9)

        available_time = params.fo_max_time - elapsed_time

        if elapsed_time > params.fo_max_time
            stop_by_time_limit = true
            println("\nFO time limit exceeded")
            break
        end

        # Compare not improved rounds for theta = 0
        if !improved && params.fo_theta == 0
            not_improved_rounds += 1
        end

        # Dynamically changes window_size and step_size        
        if !improved && params.fo_theta > 0
            not_improved_rounds += 1
            window_size += params.fo_theta
            if window_size - step_size > 2
                step_size += 1
            end
        else # Reset window_size and step_size
            window_size = Int64(max(2, floor(data.NT / params.fo_window_size_factor)))
            step_size = Int64(ceil(window_size / params.fo_step_size_factor))
        end
    end

    for i in I
        for t in T
            
            set_upper_bound(model[:w][i, t], 1.0)
            set_lower_bound(model[:w][i, t], 0.0)
            set_upper_bound(model[:y][i, t], 1.0)
            set_lower_bound(model[:y][i, t], 0.0)
          
        end
    end

    # Check if solution is feasible (integer)
    integerSol = true

    for i in I, t in T
        if value(sol.y[i, t]) > 0.01 && value(sol.y[i, t]) < 0.99
            integerSol = false
            println("FFAAIILL!!")
            break
        end
    end

    if integerSol == true
        if stop_by_time_limit == false
            stats.sol_status = "Success"
            println("\nRF found an integer solution with value: ", sol.primal_bound)
        else
            stats.sol_status = "Success_TL"
            println("\nRF stopped by time limit but found an integer solution with value: ", sol.primal_bound)
        end
    else
        if stop_by_time_limit == false
            stats.sol_status = "Fail: !Int"
            println("\nRF could not find an integer solution")
        else
            stats.sol_status = "Fail: !Int_TL"
            println("\nRF stopped by time limit and could not find an integer solution")
        end
    end

    time_end = time_ns()

    stats.fo_time = ((time_end - time_start) * 1e-9)
    stats.fo_UB = best_UB
    stats.fo_improvement = 100 * ((stats.rf_UB - stats.fo_UB) / stats.rf_UB)

    println("\nBest primal bound = ", best_UB)

    println("\nfo_time = ", stats.fo_time)

    println("===================== Finished FO! ======================")

end
function fix_and_optimize_per_period_1!(data::InstanceData, model::Model,
    sol::Union{StdFormModelSolution, FlFormModelSolution, SrFormModelSolution}, params::ExperimentParameters, stats::StatisticsData)

    println("===================== Running FO... =====================")

    I = 1:data.NI
    T = 1:data.NT
    

    elapsed_time = 0.0

    # Get fix-and-optimize start time
    time_start = time_ns()

    stats.fo_iterations = 0
    stats.fo_improvement = 0.0
    stats.fo_time = 0.0

    stop_by_time_limit = false

    best_UB = sol.primal_bound
    println("Current primal bound = ", best_UB)
    curr_UB = best_UB
    improved = true
    improvement = 0.0

    # Inititalize window
    window_size = Int64(max(2, round(data.NT / params.fo_window_size_factor)))
    step_size = Int64(round(window_size / params.fo_step_size_factor))
    max_window_size = Int64(round(data.NT / params.fo_max_window_size_factor))#  params.fo_max_window_size_factor * data.NT




    if params.fo_theta == 0
        params.fo_max_rounds_without_improvement = 1
    end

    # sol_is_feasible = Formulation.check_solution(data, sol)

    # if sol_is_feasible == true
    #     println("Solution is FEASIBLE")
    #     for k in M, t in T
    #         set_upper_bound(model[:cap_viol_var][t], 0.0)
    #     end
    # else
    #     println("Solution is INFEASIBLE")
    # end

    not_improved_rounds = 0
    stats.fo_rounds = 0
    params.fo_max_time = params.total_time_limit - stats.rf_time
    available_time = params.fo_max_time  
    println("available time: ", available_time, " RF solution = ", best_UB)

    Parameters.set_fo_parameters(model, params)

    while elapsed_time < params.fo_max_time && window_size <= max_window_size && not_improved_rounds < params.fo_max_rounds_without_improvement
        # println("rodar a primeira parte")
        stats.fo_rounds += 1

        # Inititalize window
        window_begin = 1

        # Compute end of window
        window_end = window_begin + window_size - 1

        improved = false
        improvement = 0.0

        #=  for t in window_begin:window_end
            for i in I
                println("Upper bound of w[$i, $t] = ", upper_bound(model[:w][i, t]))
            end
        end =#

        while window_end <= data.NT && window_begin <= data.NT && elapsed_time < params.fo_max_time  && sol.status != INFEASIBLE

            stats.fo_iterations += 1
            print( "\nrnd = ", stats.fo_rounds, " iter: ",  stats.fo_iterations,  " [$window_begin..$window_end]", " best_UB = ", round(best_UB, digits = 5))

            # Fix variables before the window
            if window_begin > 1
                for t in 1:window_begin - 1
                    for i in I
                        
                        if sol.y[i, t] >= 0.9
                            set_lower_bound(model[:y][i, t], 1.0)
                            set_upper_bound(model[:y][i, t], 1.0)
                        end
                        if sol.y[i, t] <= 0.1
                            set_lower_bound(model[:y][i, t], 0.0)
                            set_upper_bound(model[:y][i, t], 0.0)
                        end
                        if sol.w[i, t] >= 0.9
                            set_lower_bound(model[:w][i, t], 1.0)
                            set_upper_bound(model[:w][i, t], 1.0)
                        end
                        if sol.w[i, t] <= 0.1
                            set_lower_bound(model[:w][i, t], 0.0)
                            set_upper_bound(model[:w][i, t], 0.0)
                        end
                        
                    end
                end
            end

            # Define variables as binary inside the window
            for t in window_begin:window_end
                for i in I
                    
                    if !is_binary(model[:y][i, t])
                        set_binary(model[:y][i, t])
                    end
                    set_lower_bound(model[:y][i, t], 0.0)
                    set_upper_bound(model[:y][i, t], 1.0)

                    if !is_binary(model[:w][i, t])
                        set_binary(model[:w][i, t])
                    end
                    set_lower_bound(model[:w][i, t], 0.0)
                    set_upper_bound(model[:w][i, t], 1.0)

                    # set_start_value(model[:y][i,k,t], sol.y[i,k,t])
                    
                end
            end

            # Fix variables after the window
            if window_begin == 1
                for t in window_end + 1:data.NT
                    for i in I

                        if sol.y[i, t] >= 0.9
                            set_lower_bound(model[:y][i, t], 1.0)
                            set_upper_bound(model[:y][i, t], 1.0)
                        end
                        if sol.y[i, t] <= 0.1
                            set_lower_bound(model[:y][i, t], 0.0)
                            set_upper_bound(model[:y][i, t], 0.0)
                        end
                        if sol.w[i, t] >= 0.9
                            set_lower_bound(model[:w][i, t], 1.0)
                            set_upper_bound(model[:w][i, t], 1.0)
                        end
                        if sol.w[i, t] <= 0.1
                            set_lower_bound(model[:w][i, t], 0.0)
                            set_upper_bound(model[:w][i, t], 0.0)
                        end

                    end
                end
            end

            # Solve restricted model
            set_time_limit_sec(model, min(available_time, params.fo_restricted_model_max_time))
            #write_to_file(model, "model.lp")
            set_optimizer_attributes(model, "NodeLimit" => 0)

            if params.MIP_model == "std_form"
                println("Solving std_form model")
                Formulation.solve_stdform_model!(model, data, sol, stats)
            elseif params.MIP_model == "fac_loc_form"
                println("Solving fl_form model")
                ReformulationFL.solve_flform_model!(model, data, sol, stats)
            elseif params.MIP_model == "short_path_form"
                println("Solving rs_form model")
                ReformulationSR.solve_rsform_model!(model, data, sol, stats)
            end

            if termination_status(model) != MOI.INFEASIBLE  
                if !has_values(model)
                    set_optimizer_attributes(model, "NodeLimit" => Inf )
                    
                    set_optimizer_attributes(model, "SolutionLimit" => 1)
    
                    if params.MIP_model == "std_form"
                        println("Solving std_form model")
                        Formulation.solve_stdform_model!(model, data, sol, stats)
                    elseif params.MIP_model == "fac_loc_form"
                        println("Solving fl_form model")
                        ReformulationFL.solve_flform_model!(model, data, sol, stats)
                    elseif params.MIP_model == "short_path_form"
                        println("Solving rs_form model")
                        ReformulationSR.solve_rsform_model!(model, data, sol, stats)
                    end
    
                    set_optimizer_attributes(model, "SolutionLimit" =>    2000000000)
                end
                
            end
             
            sol.status = termination_status(model)
            println(sol.status)
            curr_UB = sol.primal_bound
            improvement = 0.0

            if best_UB - curr_UB > 0.000001
                improved = true
                not_improved_rounds = 0
                improvement = 100 * ((best_UB - curr_UB) / best_UB)
                best_UB = curr_UB
            end

            print(" curr_obj = ", round(curr_UB, digits = 5),  " improvement = ", round(improvement, digits = 4))

            curr_time = time_ns()
            elapsed_time = ((curr_time - time_start) * 1e-9)
            print(" elapsed_time = ", round(elapsed_time, digits = 2))

            available_time = params.fo_max_time - elapsed_time
            print(" available time: ", round(available_time, digits = 4))

            if elapsed_time > params.fo_max_time
                stop_by_time_limit = true
                println("\nFO time limit exceeded")
                break
            end
            
            if window_end == data.NT
                break
            end

            # Update FO window
            window_begin += step_size
            window_end = min(window_begin + window_size - 1, data.NT)

        end

        println("\nSolution was improved in round $(stats.fo_rounds)? R.: ", improved)
        curr_time = time_ns()
        elapsed_time = ((curr_time - time_start) * 1e-9)

        available_time = params.fo_max_time - elapsed_time

        if elapsed_time > params.fo_max_time
            stop_by_time_limit = true
            println("\nFO time limit exceeded")
            break
        end

        # Compare not improved rounds for theta = 0
        if !improved && params.fo_theta == 0
            not_improved_rounds += 1
        end

        # Dynamically changes window_size and step_size        
        if !improved && params.fo_theta > 0
            not_improved_rounds += 1
            window_size += params.fo_theta
            if window_size - step_size > 2
                step_size += 1
            end
        else # Reset window_size and step_size
            window_size = Int64(max(2, floor(data.NT / params.fo_window_size_factor)))
            step_size = Int64(ceil(window_size / params.fo_step_size_factor))
        end
    end
    set_optimizer_attributes(model, "NodeLimit" => Inf )
    # Check if solution is feasible (integer)
    integerSol = true

    for i in I, t in T
        if value(sol.y[i, t]) > 0.01 && value(sol.y[i, t]) < 0.99
            integerSol = false
            println("FFAAIILL!!")
            break
        end
    end

    if integerSol == true
        if stop_by_time_limit == false
            stats.sol_status = "Success"
            println("\nRF found an integer solution with value: ", sol.primal_bound)
        else
            stats.sol_status = "Success_TL"
            println("\nRF stopped by time limit but found an integer solution with value: ", sol.primal_bound)
        end
    else
        if stop_by_time_limit == false
            stats.sol_status = "Fail: !Int"
            println("\nRF could not find an integer solution")
        else
            stats.sol_status = "Fail: !Int_TL"
            println("\nRF stopped by time limit and could not find an integer solution")
        end
    end


    for i in I
        for t in T
            
            set_upper_bound(model[:w][i, t], 1.0)
            set_lower_bound(model[:w][i, t], 0.0)
            set_upper_bound(model[:y][i, t], 1.0)
            set_lower_bound(model[:y][i, t], 0.0)
          
        end
    end
   # write_to_file(model, "model1.lp")
    time_end = time_ns()

    stats.fo_time = ((time_end - time_start) * 1e-9)
    stats.fo_UB = best_UB
    stats.fo_improvement = 100 * ((stats.rf_UB - stats.fo_UB) / stats.rf_UB)

    println("\nBest primal bound = ", best_UB)

    println("\nfo_time = ", stats.fo_time)

    println("===================== Finished FO! ======================")

end
function fix_and_optimize_per_period_Relaxation_Search!(data::InstanceData, model::Model,
    sol::Union{StdFormModelSolution, FlFormModelSolution, SrFormModelSolution}, params::ExperimentParameters, stats::StatisticsData)

    println("===================== Running FO... =====================")

    I = 1:data.NI
    T = 1:data.NT
    

    elapsed_time = 0.0

    # Get fix-and-optimize start time
    time_start = time_ns()

    stats.fo_iterations = 0
    stats.fo_improvement = 0.0
    stats.fo_time = 0.0

    stop_by_time_limit = false

    best_UB = sol.primal_bound
    println("Current primal bound = ", best_UB)
    curr_UB = best_UB
    improved = true
    improvement = 0.0

    # Inititalize window
    window_size = Int64(max(2, round(data.NT / params.fo_window_size_factor)))
    step_size = Int64(round(window_size / params.fo_step_size_factor))
    max_window_size = Int64(round(data.NT / params.fo_max_window_size_factor))#  params.fo_max_window_size_factor * data.NT

    


    #=  if params.fo_theta == 0
        params.fo_max_rounds_without_improvement = 1
    end =#

    # sol_is_feasible = Formulation.check_solution(data, sol)

    # if sol_is_feasible == true
    #     println("Solution is FEASIBLE")
    #     for k in M, t in T
    #         set_upper_bound(model[:cap_viol_var][t], 0.0)
    #     end
    # else
    #     println("Solution is INFEASIBLE")
    # end

    not_improved_rounds = 0
    stats.fo_rounds = 0
    params.fo_max_time = params.total_time_limit - stats.rf_time
    available_time = params.fo_max_time  
    println("available time: ", available_time, " RF solution = ", best_UB)

    Parameters.set_fo_parameters(model, params)

    # println("rodar a primeira parte")
    stats.fo_rounds += 1

    # Inititalize window
    window_begin = 1

    # Compute end of window
    window_end = window_begin + window_size - 1

    improved = false
    improvement = 0.0
    solutionsrelex = []
    solutionsrelexwb = []
    solutionsrelexwe = []
    y_copy = Dict((i, t) => sol.y[i, t] for i in I, t in T)
    w_copy = Dict((i, t) => sol.w[i, t] for i in I, t in T)
    window_antt = 0
    window_plus = 0

    while elapsed_time < params.fo_max_time  && not_improved_rounds <  params.fo_max_rounds_without_improvement

        
       
        

        if not_improved_rounds == 0
            println("relaxation search")
            # Inititalize window
            window_begin = 1

            # Compute end of window
            window_end = window_begin + window_size - 1

            solutionsrelex = []
            solutionsrelexwb = []
            solutionsrelexwe = []

            y_copy = Dict((i, t) => sol.y[i, t] for i in I, t in T)
            w_copy = Dict((i, t) => sol.w[i, t] for i in I, t in T)
            
           #=  for t in T
                for i in I
                  println("y_copy[$i, $t] = ", y_copy[i,t])
                end
            end =#

            while window_begin <= window_end - step_size
                
               # println("window_begin = ", window_begin)
                #println("window_end = ", window_end)
             # println("window_size = ", window_size, " step_size = ", step_size, " max_window_size = ", max_window_size)
                # Restore model from model_copy
                
                    for i in I
                        for t in T
                            if !is_binary(model[:y][i,t])
                                set_binary(model[:y][i,t])
                                set_binary(model[:w][i,t])
                            end
                            if !is_fixed(model[:y][i,t])
                                fix(model[:y][i,t], y_copy[i,t], force = true)
                                fix(model[:w][i,t], w_copy[i,t], force = true)
                            end
                    
                        end
                    end
                
            
                # Define variables as continuas inside the window
                for t in window_begin:window_end
                    for i in I
                        
                        if is_fixed(model[:y][i,t])  
                            
                            
                            unfix(model[:y][i,t])
                            unfix(model[:w][i,t])
                        end

                    
                        if is_binary(model[:y][i,t])                      
                                
                            unset_binary(model[:y][i,t])
                            set_lower_bound(model[:y][i,t], 0.0)
                            set_upper_bound(model[:y][i,t], 1.0)

                            unset_binary(model[:w][i,t])
                            set_lower_bound(model[:w][i,t], 0.0)
                            set_upper_bound(model[:w][i,t], 1.0)
                            
                        end   
                    end
                end


                if params.MIP_model == "std_form"
                   # println("Solving std_form model")
                    Formulation.solve_stdform_model!(model, data, sol, stats)
                elseif params.MIP_model == "fac_loc_form"
                    println("Solving fl_form model")
                    ReformulationFL.solve_flform_model!(model, data, sol, stats)
                elseif params.MIP_model == "short_path_form"
                    println("Solving rs_form model")
                    ReformulationSR.solve_rsform_model!(model, data, sol, stats)
                end
                if has_values(model)
                    soluRelex = objective_value(model)
                else
                    soluRelex = Inf
                end
                
            
             # println("solução = ", soluRelex)
                push!(solutionsrelex, soluRelex)
                push!(solutionsrelexwb, window_begin)
                push!(solutionsrelexwe, window_end)

                # Update FO window
                window_begin += step_size
                window_end = min(window_begin +  window_size - 1 , data.NT)
                
                
            end


           
            window_plus = 0




        end


        window_diff = false
        
        while window_diff == false && window_begin <= data.NT - step_size
            println("calcular a janela")
           

            # println("Solutions: ", solutionsrelex)
            solutionsrelexord = sort(solutionsrelex)
            #println("Solutions ordered: ", solutionsrelexord)
            possicao = not_improved_rounds + 1 + window_plus
            println("Position: ", possicao)
            min_solution = solutionsrelexord[possicao]
            
            index_min_solution = findfirst(==(min_solution), solutionsrelex)
            #  println("Minimum solution value: ", min_solution)
            #index_min_solution = argmin(solutionsrelex)
            # println("Index of minimum solution: ", index_min_solution)

            window_begin = solutionsrelexwb[index_min_solution]
            window_end = solutionsrelexwe[index_min_solution]

            if window_begin == window_antt
                println("window_begin = window_antt")
                window_plus += 1
            else
                window_diff = true
            end
        end
        window_antt = window_begin
        println("Window begin: ", window_begin)
        println("Window end: ", window_end)

        
        #=  println("ycopy antes do FO= ", y_copy)
        println("wcopy antes do FO= ", w_copy) =#

        for i in I
            for t in T
                if !is_binary(model[:y][i,t])
                    set_binary(model[:y][i,t])
                    set_binary(model[:w][i,t])
                end
                #if !is_fixed(model[:y][i,t])
                    fix(model[:y][i,t], y_copy[i,t], force = true)
                    fix(model[:w][i,t], w_copy[i,t], force = true)
               # end
           
            end
        end
        # Define variables as intege inside the window
        for t in window_begin:window_end
            for i in I
                
                if is_fixed(model[:y][i,t])  
                    unfix(model[:y][i,t])
                    unfix(model[:w][i,t])
                end
            
                 

                set_lower_bound(model[:y][i,t], 0.0)
                set_upper_bound(model[:y][i,t], 1.0) 
                set_lower_bound(model[:w][i,t], 0.0)
                set_upper_bound(model[:w][i,t], 1.0)
            end
        end



        #=  # check if there are no fixed or continuous values left
        for t in T
            for i in I
                if  is_fixed(model[:y][i,t])
    
                    unfix(model[:y][i,t])
                    unfix(model[:w][i,t])
                end
              
                if !is_binary(model[:y][i,t])
                    set_binary(model[:y][i,t])
                    set_binary(model[:w][i,t])
                end
            end
        end =#

        #= # Fix variables before the window
        if window_begin > 1
            for t in 1:window_begin - 1
                for i in I
                    
                    if sol.y[i, t] >= 0.9
                        set_lower_bound(model[:y][i, t], 1.0)
                        set_upper_bound(model[:y][i, t], 1.0)
                    end
                    if sol.y[i, t] <= 0.1
                        set_lower_bound(model[:y][i, t], 0.0)
                        set_upper_bound(model[:y][i, t], 0.0)
                    end
                    if sol.w[i, t] >= 0.9
                        set_lower_bound(model[:w][i, t], 1.0)
                        set_upper_bound(model[:w][i, t], 1.0)
                    end
                    if sol.w[i, t] <= 0.1
                        set_lower_bound(model[:w][i, t], 0.0)
                        set_upper_bound(model[:w][i, t], 0.0)
                    end
                    
                end
            end
        end

        # Define variables as binary inside the window
        for t in window_begin:window_end
            for i in I
                
                if !is_binary(model[:y][i, t])
                    set_binary(model[:y][i, t])
                end
                set_lower_bound(model[:y][i, t], 0.0)
                set_upper_bound(model[:y][i, t], 1.0)

                if !is_binary(model[:w][i, t])
                    set_binary(model[:w][i, t])
                end
                set_lower_bound(model[:w][i, t], 0.0)
                set_upper_bound(model[:w][i, t], 1.0)

                
            end
        end

        # Fix variables after the window
        if window_begin == 1
            for t in window_end + 1:data.NT
                for i in I

                    if sol.y[i, t] >= 0.9
                        set_lower_bound(model[:y][i, t], 1.0)
                        set_upper_bound(model[:y][i, t], 1.0)
                    end
                    if sol.y[i, t] <= 0.1
                        set_lower_bound(model[:y][i, t], 0.0)
                        set_upper_bound(model[:y][i, t], 0.0)
                    end
                    if sol.w[i, t] >= 0.9
                        set_lower_bound(model[:w][i, t], 1.0)
                        set_upper_bound(model[:w][i, t], 1.0)
                    end
                    if sol.w[i, t] <= 0.1
                        set_lower_bound(model[:w][i, t], 0.0)
                        set_upper_bound(model[:w][i, t], 0.0)
                    end

                end
            end
        end =#

        
        # Solve restricted model
        set_time_limit_sec(model, min(available_time, params.fo_restricted_model_max_time))

        # write_to_file(model, "model.lp")

        if params.MIP_model == "std_form"
            println("Solving std_form model")
            Formulation.solve_stdform_model!(model, data, sol, stats)
        elseif params.MIP_model == "fac_loc_form"
            println("Solving fl_form model")
            ReformulationFL.solve_flform_model!(model, data, sol, stats)
        elseif params.MIP_model == "short_path_form"
            println("Solving rs_form model")
            ReformulationSR.solve_rsform_model!(model, data, sol, stats)
        end
            
        sol.status = termination_status(model)
        println(sol.status)
        curr_UB = sol.primal_bound
        println("Current primal bound = ", curr_UB)
        improvement = 0.0

        if best_UB - curr_UB > 0.01
            improved = true
            not_improved_rounds = 0
            improvement = 100 * ((best_UB - curr_UB) / best_UB)
            best_UB = curr_UB
        else
            not_improved_rounds += 1
        end
        println("not_improved_rounds = ", not_improved_rounds)

        print(" curr_obj = ", round(curr_UB, digits = 5),  " improvement = ", round(improvement, digits = 4))

        curr_time = time_ns()
        elapsed_time = ((curr_time - time_start) * 1e-9)
        print(" elapsed_time = ", round(elapsed_time, digits = 2))

        available_time = params.fo_max_time - elapsed_time
        print(" available time: ", round(available_time, digits = 4))

        if elapsed_time > params.fo_max_time
            stop_by_time_limit = true
            println("\nFO time limit exceeded")
            break
        end
        
    end


        #= 
        while window_end <= data.NT && window_begin <= data.NT && elapsed_time < params.fo_max_time  && sol.status != INFEASIBLE

            stats.fo_iterations += 1
            print( "\nrnd = ", stats.fo_rounds, " iter: ",  stats.fo_iterations,  " [$window_begin..$window_end]", " best_UB = ", round(best_UB, digits = 5))

           

            

            

        end

        println("\nSolution was improved in round $(stats.fo_rounds)? R.: ", improved)
        curr_time = time_ns()
        elapsed_time = ((curr_time - time_start) * 1e-9)

        available_time = params.fo_max_time - elapsed_time

        if elapsed_time > params.fo_max_time
            stop_by_time_limit = true
            println("\nFO time limit exceeded")
            break
        end

        # Compare not improved rounds for theta = 0
        if !improved && params.fo_theta == 0
            not_improved_rounds += 1
        end

        # Dynamically changes window_size and step_size        
        if !improved && params.fo_theta > 0
            not_improved_rounds += 1
            window_size += params.fo_theta
            if window_size - step_size > 2
                step_size += 1
            end
        else # Reset window_size and step_size
            window_size = Int64(max(2, floor(data.NT / params.fo_window_size_factor)))
            step_size = Int64(ceil(window_size / params.fo_step_size_factor))
        end
       
    end =#
  
    # Check if solution is feasible (integer)
    integerSol = true

    for i in I, t in T
        if value(sol.y[i, t]) > 0.01 && value(sol.y[i, t]) < 0.99
            integerSol = false
            println("FFAAIILL!!")
            break
        end
    end

    if integerSol == true
        if stop_by_time_limit == false
            stats.sol_status = "Success"
            println("\nRF found an integer solution with value: ", sol.primal_bound)
        else
            stats.sol_status = "Success_TL"
            println("\nRF stopped by time limit but found an integer solution with value: ", sol.primal_bound)
        end
    else
        if stop_by_time_limit == false
            stats.sol_status = "Fail: !Int"
            println("\nRF could not find an integer solution")
        else
            stats.sol_status = "Fail: !Int_TL"
            println("\nRF stopped by time limit and could not find an integer solution")
        end
    end

    time_end = time_ns()

    stats.fo_time = ((time_end - time_start) * 1e-9)
    stats.fo_UB = best_UB
    stats.fo_improvement = 100 * ((stats.rf_UB - stats.fo_UB) / stats.rf_UB)

    println("\nBest primal bound = ", best_UB)

    println("\nfo_time = ", stats.fo_time)

 
    println("===================== Finished FO! ======================")

end
function fix_and_optimize_per_period_Relaxation_Search2!(data::InstanceData, model::Model,
    sol::Union{StdFormModelSolution, FlFormModelSolution, SrFormModelSolution}, params::ExperimentParameters, stats::StatisticsData)

    println("===================== Running FO... =====================")

    I = 1:data.NI
    T = 1:data.NT
    

    elapsed_time = 0.0

    # Get fix-and-optimize start time
    time_start = time_ns()

    stats.fo_iterations = 0
    stats.fo_improvement = 0.0
    stats.fo_time = 0.0

    stop_by_time_limit = false

    best_UB = sol.primal_bound
    println("Current primal bound = ", best_UB)
    curr_UB = best_UB
    improved = true
    improvement = 0.0

    # Inititalize window
    window_size = Int64(max(2, round(data.NT / params.fo_window_size_factor)))
    step_size = Int64(round(window_size / params.fo_step_size_factor))
    max_window_size = Int64(round(data.NT / params.fo_max_window_size_factor))#  params.fo_max_window_size_factor * data.NT

    


    if params.fo_theta == 0
        params.fo_max_rounds_without_improvement = 1
    end

    # sol_is_feasible = Formulation.check_solution(data, sol)

    # if sol_is_feasible == true
    #     println("Solution is FEASIBLE")
    #     for k in M, t in T
    #         set_upper_bound(model[:cap_viol_var][t], 0.0)
    #     end
    # else
    #     println("Solution is INFEASIBLE")
    # end

    not_improved_rounds = 0
    stats.fo_rounds = 0
    params.fo_max_time = params.total_time_limit - stats.rf_time
    available_time = params.fo_max_time  
    println("available time: ", available_time, " RF solution = ", best_UB)

    Parameters.set_fo_parameters(model, params)

    # println("rodar a primeira parte")
    stats.fo_rounds += 1

    # Inititalize window
    window_begin = 1

    # Compute end of window
    window_end = window_begin + window_size - 1

    #improved = false
    improvement = 0.0
    solutionsrelex = []
    solutionsrelexwb = []
    solutionsrelexwe = []
    #= y_copy = Dict((i, t) => sol.y[i, t] for i in I, t in T)
    w_copy = Dict((i, t) => sol.w[i, t] for i in I, t in T) =#
    

    while elapsed_time < params.fo_max_time  && improved == true

        
            improved = false
        

        
            println("rodar relaxações")
            # Inititalize window
            window_begin = 1

            # Compute end of window
            window_end = window_begin + window_size - 1

            solutionsrelex = []
            solutionsrelexwb = []
            solutionsrelexwe = []

            y_copy = Dict((i, t) => sol.y[i, t] for i in I, t in T)
            w_copy = Dict((i, t) => sol.w[i, t] for i in I, t in T)

            #= for i in I
                for t in T
                  println("y_copy[$i, $t] = ", y_copy[i,t])
                end
            end =#

            while window_begin <= window_end - step_size
                
                # println("window_begin = ", window_begin)
                #println("window_end = ", window_end)
                # println("window_size = ", window_size, " step_size = ", step_size, " max_window_size = ", max_window_size)
                # Restore model from model_copy
                
                    for i in I
                        for t in T
                            if !is_binary(model[:y][i,t])
                                set_binary(model[:y][i,t])
                                set_binary(model[:w][i,t])
                            end
                            if !is_fixed(model[:y][i,t])
                                fix(model[:y][i,t], y_copy[i,t], force = true)
                                fix(model[:w][i,t], w_copy[i,t], force = true)
                            end
                    
                        end
                    end
                
            
                # Define variables as continuas inside the window
                for t in window_begin:window_end
                    for i in I
                        
                        if is_fixed(model[:y][i,t])  
                            
                            
                            unfix(model[:y][i,t])
                            unfix(model[:w][i,t])
                        end

                    
                        if is_binary(model[:y][i,t])                      
                                
                            unset_binary(model[:y][i,t])
                            set_lower_bound(model[:y][i,t], 0.0)
                            set_upper_bound(model[:y][i,t], 1.0)

                            unset_binary(model[:w][i,t])
                            set_lower_bound(model[:w][i,t], 0.0)
                            set_upper_bound(model[:w][i,t], 1.0)
                            
                        end   
                    end
                end


                if params.MIP_model == "std_form"
                   # println("Solving std_form model")
                    Formulation.solve_stdform_model!(model, data, sol, stats)
                elseif params.MIP_model == "fac_loc_form"
                    println("Solving fl_form model")
                    ReformulationFL.solve_flform_model!(model, data, sol, stats)
                elseif params.MIP_model == "short_path_form"
                    println("Solving rs_form model")
                    ReformulationSR.solve_rsform_model!(model, data, sol, stats)
                end
                if has_values(model)
                    soluRelex = objective_value(model)
                else
                    soluRelex = Inf
                end
                
            
                # println("solução = ", soluRelex)
                push!(solutionsrelex, soluRelex)
                push!(solutionsrelexwb, window_begin)
                push!(solutionsrelexwe, window_end)

                # Update FO window
                window_begin += step_size
                window_end = min(window_begin +  window_size - 1 , data.NT)
                
                
            end
        

        # println("Solutions: ", solutionsrelex)
        solutionsrelexord = sort(solutionsrelex)
        println("Solutions ordered: ", solutionsrelexord)#=  
        possicao = not_improved_rounds +1
        min_solution = solutionsrelexord[possicao] =# 


        #=    for i in 1:length(solutionsrelexord)
            k = findfirst(==(solutionsrelexord[i]), solutionsrelex)
            window_begin = solutionsrelexwb[k]
            window_end = solutionsrelexwe[k]
            #println("For solution $i: window_begin = $window_begin, window_end = $window_end")
        end =#


        #=   index_min_solution = findfirst(==(min_solution), solutionsrelex)
        #  println("Minimum solution value: ", min_solution)
        #index_min_solution = argmin(solutionsrelex)
        # println("Index of minimum solution: ", index_min_solution)

        window_begin = solutionsrelexwb[index_min_solution]
        window_end = solutionsrelexwe[index_min_solution]
        println("Window begin: ", window_begin)
        println("Window end: ", window_end) =#

        #=  println("ycopy antes do FO= ", y_copy)
        println("wcopy antes do FO= ", w_copy) =#

        println("rodar principal")


        for k in 1:length(solutionsrelexord)

            f = findfirst(==(solutionsrelexord[k]), solutionsrelex)
            window_begin = solutionsrelexwb[f]
            window_end = solutionsrelexwe[f]

            println("Window begin: ", window_begin)
            println("Window end: ", window_end)

            for i in I
                for t in T
                    if !is_binary(model[:y][i,t])
                        set_binary(model[:y][i,t])
                        set_binary(model[:w][i,t])
                    end
                    if k == 1

                        fix(model[:y][i,t], y_copy[i,t], force = true)
                        fix(model[:w][i,t], w_copy[i,t], force = true)
                    else
                        fix(model[:y][i,t], sol.y[i,t], force = true)
                        fix(model[:w][i,t], sol.w[i,t], force = true)

                    end
            
                end
            end
            # Define variables as intege inside the window
            for t in window_begin:window_end
                for i in I
                    
                    if is_fixed(model[:y][i,t])  
                        unfix(model[:y][i,t])
                        unfix(model[:w][i,t])
                    end
                
                    

                    set_lower_bound(model[:y][i,t], 0.0)
                    set_upper_bound(model[:y][i,t], 1.0) 
                    set_lower_bound(model[:w][i,t], 0.0)
                    set_upper_bound(model[:w][i,t], 1.0)
                end
            end



            #=  # check if there are no fixed or continuous values left
            for t in T
                for i in I
                    if  is_fixed(model[:y][i,t])
        
                        unfix(model[:y][i,t])
                        unfix(model[:w][i,t])
                    end
                
                    if !is_binary(model[:y][i,t])
                        set_binary(model[:y][i,t])
                        set_binary(model[:w][i,t])
                    end
                end
            end =#

            #= # Fix variables before the window
            if window_begin > 1
                for t in 1:window_begin - 1
                    for i in I
                        
                        if sol.y[i, t] >= 0.9
                            set_lower_bound(model[:y][i, t], 1.0)
                            set_upper_bound(model[:y][i, t], 1.0)
                        end
                        if sol.y[i, t] <= 0.1
                            set_lower_bound(model[:y][i, t], 0.0)
                            set_upper_bound(model[:y][i, t], 0.0)
                        end
                        if sol.w[i, t] >= 0.9
                            set_lower_bound(model[:w][i, t], 1.0)
                            set_upper_bound(model[:w][i, t], 1.0)
                        end
                        if sol.w[i, t] <= 0.1
                            set_lower_bound(model[:w][i, t], 0.0)
                            set_upper_bound(model[:w][i, t], 0.0)
                        end
                        
                    end
                end
            end

            # Define variables as binary inside the window
            for t in window_begin:window_end
                for i in I
                    
                    if !is_binary(model[:y][i, t])
                        set_binary(model[:y][i, t])
                    end
                    set_lower_bound(model[:y][i, t], 0.0)
                    set_upper_bound(model[:y][i, t], 1.0)

                    if !is_binary(model[:w][i, t])
                        set_binary(model[:w][i, t])
                    end
                    set_lower_bound(model[:w][i, t], 0.0)
                    set_upper_bound(model[:w][i, t], 1.0)

                    
                end
            end

            # Fix variables after the window
            if window_begin == 1
                for t in window_end + 1:data.NT
                    for i in I

                        if sol.y[i, t] >= 0.9
                            set_lower_bound(model[:y][i, t], 1.0)
                            set_upper_bound(model[:y][i, t], 1.0)
                        end
                        if sol.y[i, t] <= 0.1
                            set_lower_bound(model[:y][i, t], 0.0)
                            set_upper_bound(model[:y][i, t], 0.0)
                        end
                        if sol.w[i, t] >= 0.9
                            set_lower_bound(model[:w][i, t], 1.0)
                            set_upper_bound(model[:w][i, t], 1.0)
                        end
                        if sol.w[i, t] <= 0.1
                            set_lower_bound(model[:w][i, t], 0.0)
                            set_upper_bound(model[:w][i, t], 0.0)
                        end

                    end
                end
            end =#

            
            # Solve restricted model
            set_time_limit_sec(model, min(available_time, params.fo_restricted_model_max_time))

            # write_to_file(model, "model.lp")

            if params.MIP_model == "std_form"
                println("Solving std_form model")
                Formulation.solve_stdform_model!(model, data, sol, stats)
            elseif params.MIP_model == "fac_loc_form"
                println("Solving fl_form model")
                ReformulationFL.solve_flform_model!(model, data, sol, stats)
            elseif params.MIP_model == "short_path_form"
                println("Solving rs_form model")
                ReformulationSR.solve_rsform_model!(model, data, sol, stats)
            end

            sol.status = termination_status(model)
            println(sol.status)
            curr_UB = sol.primal_bound
            improvement = 0.0

            if best_UB - curr_UB > 0.000001
                improved = true
                not_improved_rounds = 0
                improvement = 100 * ((best_UB - curr_UB) / best_UB)
                best_UB = curr_UB
            end

            print(" curr_obj = ", round(curr_UB, digits = 5),  " improvement = ", round(improvement, digits = 4))

            curr_time = time_ns()
            elapsed_time = ((curr_time - time_start) * 1e-9)
            print(" elapsed_time = ", round(elapsed_time, digits = 2))

            available_time = params.fo_max_time - elapsed_time
            print(" available time: ", round(available_time, digits = 4))

            if elapsed_time > params.fo_max_time
                stop_by_time_limit = true
                println("\nFO time limit exceeded")
                break
            end
            
            if window_end == data.NT
                break
            end

            # Update FO window
            window_begin += step_size
            window_end = min(window_begin + window_size - 1, data.NT)


        end        

       
        println("\nSolution was improved in round $(stats.fo_rounds)? R.: ", improved)
        curr_time = time_ns()
        elapsed_time = ((curr_time - time_start) * 1e-9)

        available_time = params.fo_max_time - elapsed_time

        if elapsed_time > params.fo_max_time
            stop_by_time_limit = true
            println("\nFO time limit exceeded")
            break
        end
      
        # Compare not improved rounds for theta = 0
        if !improved && params.fo_theta == 0
            not_improved_rounds += 1
        end

        # Dynamically changes window_size and step_size        
        if !improved && params.fo_theta > 0
            not_improved_rounds += 1
            window_size += params.fo_theta
            if window_size - step_size > 2
                step_size += 1
            end
        else # Reset window_size and step_size
            window_size = Int64(max(2, floor(data.NT / params.fo_window_size_factor)))
            step_size = Int64(ceil(window_size / params.fo_step_size_factor))
        end
       #=  for i in I
            for t in T
              println("sol.y[$i, $t] = ", sol.y[i,t])
            end
        end =#
    end


        #= 
        while window_end <= data.NT && window_begin <= data.NT && elapsed_time < params.fo_max_time  && sol.status != INFEASIBLE

            stats.fo_iterations += 1
            print( "\nrnd = ", stats.fo_rounds, " iter: ",  stats.fo_iterations,  " [$window_begin..$window_end]", " best_UB = ", round(best_UB, digits = 5))

           

            

            

        end

        println("\nSolution was improved in round $(stats.fo_rounds)? R.: ", improved)
        curr_time = time_ns()
        elapsed_time = ((curr_time - time_start) * 1e-9)

        available_time = params.fo_max_time - elapsed_time

        if elapsed_time > params.fo_max_time
            stop_by_time_limit = true
            println("\nFO time limit exceeded")
            break
        end

        # Compare not improved rounds for theta = 0
        if !improved && params.fo_theta == 0
            not_improved_rounds += 1
        end

        # Dynamically changes window_size and step_size        
        if !improved && params.fo_theta > 0
            not_improved_rounds += 1
            window_size += params.fo_theta
            if window_size - step_size > 2
                step_size += 1
            end
        else # Reset window_size and step_size
            window_size = Int64(max(2, floor(data.NT / params.fo_window_size_factor)))
            step_size = Int64(ceil(window_size / params.fo_step_size_factor))
        end
       
    end =#
  
    # Check if solution is feasible (integer)
    integerSol = true

    for i in I, t in T
        if value(sol.y[i, t]) > 0.01 && value(sol.y[i, t]) < 0.99
            integerSol = false
            println("FFAAIILL!!")
            break
        end
    end

    if integerSol == true
        if stop_by_time_limit == false
            stats.sol_status = "Success"
            println("\nRF found an integer solution with value: ", sol.primal_bound)
        else
            stats.sol_status = "Success_TL"
            println("\nRF stopped by time limit but found an integer solution with value: ", sol.primal_bound)
        end
    else
        if stop_by_time_limit == false
            stats.sol_status = "Fail: !Int"
            println("\nRF could not find an integer solution")
        else
            stats.sol_status = "Fail: !Int_TL"
            println("\nRF stopped by time limit and could not find an integer solution")
        end
    end

    time_end = time_ns()

    stats.fo_time = ((time_end - time_start) * 1e-9)
    stats.fo_UB = best_UB
    stats.fo_improvement = 100 * ((stats.rf_UB - stats.fo_UB) / stats.rf_UB)

    println("\nBest primal bound = ", best_UB)

    println("\nfo_time = ", stats.fo_time)

 
    println("===================== Finished FO! ======================")

end

end #module    
