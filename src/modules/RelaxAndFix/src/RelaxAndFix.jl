module RelaxAndFix

using JuMP
using Gurobi
using Data
using Formulation
using OutputStatistics


using Parameters
using CPUTime


export relax_and_fix

function relax_and_fix!(model::Model, data::InstanceData, data_type::InstanceTypeData, sol::StdFormModelSolution, stats::StatisticsData,params::ExperimentParameters)

    println("===================== Running RF... =====================")

    N = 1:data_type.DIMENSION

  

    elapsed_time = 0.0
    params.rf_max_time = min(params.rf_max_time, params.total_time_limit)

    # Get relax-and-fix start time
    time_start = time_ns()

    stats.rf_iterations = 0
    stats.rf_backtrack_attempts = 0
    stats.rf_time = 0.0

    stop_by_time_limit = false

    # Inititalize window
    window_begin = 1

 
    Number_of_Int = params.rf_Number_of_Int
    Number_of_fix =  params.rf_Number_of_fix
    Number_of_values = data_type.DIMENSION * data_type.DIMENSION
    
    
    Number_of_values_fixed = 0
    Number_of_values_int = 0
    # Compute end of window
    

    available_time = params.rf_max_time  
    println("available time: ", available_time)
    println("Number_of_Int: ", Number_of_Int)
    println("Number_of_fix: ", Number_of_fix)
    println("Number_of_values: ", Number_of_values)

    Parameters.set_rf_parameters(model, params)

    integerSol = false

    # Relax-and-fix main loop

    fix_i_x = []
    fix_j_x = []

    while (Number_of_values_int < Number_of_values && elapsed_time < params.rf_max_time) && integerSol == false

        stats.rf_iterations += 1
        println( "\niter: ",  stats.rf_iterations,  " Number_of_values_fixed = ", Number_of_values_fixed)

        if stats.rf_iterations == 1

            relax_integrality(model)

            println("Solving std_form model")
            Formulation.solve_stdform_model!(model, data, data_type, sol, stats)    

            
            sol.status = termination_status(model)
            
            if sol.status == MOI.INFEASIBLE
                println("\nrelaxation infeasible")
                stats.sol_status = "Fail: !Int_TL"
                break
            end

            for i in N
                for j in N
                    if i == j
                        set_binary(model[:x][i,j])
                        set_lower_bound(model[:x][i, j], 0.0)
                        set_upper_bound(model[:x][i, j], 0.0)
                        Number_of_values_int += 1
                    end
                end
            end
      #=       for i in 1:data_type.DIMENSION


                        if !is_binary(model[:u][i])
                            set_binary(model[:u][i])
                            set_lower_bound(model[:u][i], 0.0)
                        end

            end =#
        end

        

        #fix binary
        if stats.rf_iterations > 1                      
            
            N_fix = 0
            
            for i in N
                for j in N
                    if i != j
                        if is_binary(model[:x][i,j])
                            if lower_bound(model[:x][i, j]) != upper_bound(model[:x][i, j])
                                if  N_fix  < Number_of_fix

                                    if value(sol.x[i,j]) >= 0.9
                                        set_lower_bound(model[:x][i, j], 1.0)
                                        set_upper_bound(model[:x][i, j], 1.0)
                                    end
                                    if value(sol.x[i,j]) <= 0.1
                                        set_lower_bound(model[:x][i, j], 0.0)
                                        set_upper_bound(model[:x][i, j], 0.0)
                                    end


                                    N_fix += 1
                                    Number_of_values_fixed += 1
                                    
                                    push!(fix_i_x, i)
                                    push!(fix_j_x, j)                               
                                
                                end
                            end
                        end
                    end
                end
            end
            for i in N
                for j in N
                    if i != j
                        if lower_bound(model[:x][i, j]) == 1.0
                            for k in N
                                if k != i && k != j
                                    if !is_binary(model[:x][i, k])
                                        set_binary(model[:x][i, k])
                                    end
                                    set_lower_bound(model[:x][i, k], 0.0)
                                    set_upper_bound(model[:x][i, k], 0.0)
                                    Number_of_values_int += 1
                                    Number_of_values_fixed += 1

                                end
                            end


                        end
                    end
                end

            end
        end

        # Find the Number_of_Int values closest to 0 or 1 to turn into binary
        # tem de ser feito de forma diferente tem de olhar todo os valores e pegar os primeiros mais proximos de 0.5 e subtraindo 
        value_to_int = 0.005
        N_Int = 0
        
        
        if  Number_of_values -  Number_of_values_int  <  Number_of_Int
            
            Number_of_Int = Number_of_values -  Number_of_values_int  
            
        end
        println("Number_of_Int: ", Number_of_Int)
        while  N_Int  < Number_of_Int 
            for i in N
                for j in N
                    if i != j
                        if  N_Int  < Number_of_Int
                            if !is_binary(model[:x][i,j])
                                if abs(sol.x[i,j] - 1) < value_to_int 

                                    if !is_binary(model[:x][i,j])
                                        set_binary(model[:x][i,j])
                                    end
                                    set_lower_bound(model[:x][i,j], 0.0)
                                    set_upper_bound(model[:x][i,j], 1.0)
                                    
                                    N_Int += 1
                                    Number_of_values_int += 1
                                    push!(fix_i_x, i)  
                                    push!(fix_j_x, j)

                                end
                            end
                        end
                    end
                   
                end
            end

            value_to_int += 0.005

        end
        valor_1 = 0.0
        for i in N
            for j in N
                if i != j
                    if is_binary(model[:x][i,j])
                        if lower_bound(model[:x][i,j]) == 1.0
                            valor_1 += 1.0
                        end
                    end
                end
            end
        end
        println("valor_1: ", valor_1)
        #= cont_y = 0
        for j in N 
            for i in N
                if !is_binary(model[:y][i,t])
                    cont_y += 1
                end

            end

        end
        println("Number_of_Int depois: ", cont_y) =#
       
        # Solve restricted model
        set_time_limit_sec(model, min(available_time, params.rf_restricted_model_max_time))
        set_optimizer_attributes(model, "SolutionLimit" => 2000000000)
      
        
        println("Solving std_form model")
        Formulation.solve_stdform_model!(model, data, data_type, sol, stats)    

        sol.status = termination_status(model)
        #refazer para primeiro verificar se tem valores e depois verificar se é otimo ou time limit






        if sol.status == MOI.OPTIMAL || (sol.status == MOI.TIME_LIMIT && has_values(model))

            
            println(" obj = ", round(objective_value(model), digits = 2),  " ", sol.status)

            curr_time = time_ns()
            elapsed_time = ((curr_time - time_start) * 1e-9)
            println(" elapsed_time = ", round(elapsed_time, digits = 2))

            if elapsed_time > params.rf_max_time
                stop_by_time_limit = true
                println("\nRF time limit exceeded")
                break
            end

            available_time = params.rf_max_time - elapsed_time
            println(" available time: ", round(available_time, digits = 4))
         #= else




            if has_values(model)
                sol.primal_bound = objective_value(model)
            end
            
            stats.rf_backtrack_attempts += 1

            println("backtrack start")


            backtrack = 0
           # set_time_limit_sec(model, min(available_time, 10))
           
            #torna continua as variaveis binarias 
            
            while !has_values(model) && elapsed_time < params.rf_max_time
                

                    bt_i =  backtrack  * Number_of_fix     
                    bt_f = bt_i + Number_of_fix - 1
                    
                    for p = bt_i: bt_f
                        k = Number_of_values_fixed - p

                        if lower_bound(model[:y][fix_i_x[k],fix_j_x[k]]) == upper_bound(model[:y][fix_i_x[k],fix_j_x[k]])


                            

                            if !is_binary(model[:y][fix_i_x[k],fix_j_x[k]])
                                set_binary(model[:y][fix_i_x[k],fix_j_x[k]])
                            end

                            set_lower_bound(model[:y][fix_i_x[k],fix_j_x[k]], 0.0)
                            set_upper_bound(model[:y][fix_i_x[k],fix_j_x[k]], 1.0)

                            

                            if !is_binary(model[:w][fix_i_x[k],fix_j_x[k]])
                                set_binary(model[:w][fix_i_x[k],fix_j_x[k]])
                            end
                                
                            set_lower_bound(model[:w][fix_i_x[k],fix_j_x[k]], 0.0)
                            set_upper_bound(model[:w][fix_i_x[k],fix_j_x[k]], 1.0)

                        end
                    end

                

                    set_time_limit_sec(model, min(available_time, 30))
                set_optimizer_attribute(model, "SolutionLimit", 1)
                stats.rf_backtrack__rounds += 1
                #write_to_file(model, "modelBT$backtrack.lp")
                # Solve restricted model
                
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
                println("sol.status",sol.status)
                curr_time = time_ns()
                elapsed_time = ((curr_time - time_start) * 1e-9)
                println(" elapsed_time = ", round(elapsed_time, digits = 2))

                if elapsed_time > params.rf_max_time
                    stop_by_time_limit = true
                    println("\nRF time limit exceeded")
                    break
                end
                
              

                available_time = params.fo_max_time - elapsed_time
                println(" available time: ", round(available_time, digits = 4))
               

                if !has_values(model) 
                    backtrack += 1
                    if backtrack * Number_of_fix >=  Number_of_values_fixed
                        stats.sol_status = "Fail: !Int_TL"
                        return
                    end
                    
                end
                #= if sol.status == MOI.TIME_LIMIT
                    BT_time += 1
                end =#
            end

           #fix binary backtrack
            
           if has_values(model)
                bt_i = 0     
                bt_f = Number_of_fix * (backtrack +1) -1
                
                for p = bt_i: bt_f

                    k = Number_of_values_fixed-p

                    
                            
                        if value(sol.y[fix_i_x[k],fix_j_x[k]]) >= 0.9
                            set_lower_bound(model[:y][fix_i_x[k],fix_j_x[k]], 1.0)
                            set_upper_bound(model[:y][fix_i_x[k],fix_j_x[k]], 1.0)
                        end
                        if value(sol.y[fix_i_x[k],fix_j_x[k]]) <= 0.1
                            set_lower_bound(model[:y][fix_i_x[k],fix_j_x[k]], 0.0)
                            set_upper_bound(model[:y][fix_i_x[k],fix_j_x[k]], 0.0)
                        end
                        if value(sol.w[fix_i_x[k],fix_j_x[k]]) >= 0.9
                            set_lower_bound(model[:w][fix_i_x[k],fix_j_x[k]], 1.0)
                            set_upper_bound(model[:w][fix_i_x[k],fix_j_x[k]], 1.0)
                        end
                        if value(sol.w[fix_i_x[k],fix_j_x[k]]) <= 0.1
                            set_lower_bound(model[:w][fix_i_x[k],fix_j_x[k]], 0.0)
                            set_upper_bound(model[:w][fix_i_x[k],fix_j_x[k]], 0.0)
                        end
                        
                    
                end 
            end =#
        end 
       
        integerSol = true
        for j in N, i in N
            if value(sol.x[i,j]) != 0.00 && value(sol.x[i,j]) != 1.00

                integerSol = false   
            
            end
        end
       #=  # Update RF window
        if Number_of_values_fixed == Number_of_values
            break
        end

        
        if integerSol == true && window_begin != (data.NT - rf_step_size)
            println("\nRF found an integer befor  fixed all varievel")
           # break
        end =#
        println("Number_of_values_int: ", Number_of_values_int)
       #=  for i in N
            println("i: ", i)
            for j in N
                println("x[$i,$j] = ", round(sol.x[i,j], digits = 2), " ")
            end
            
        end =#
    end
    


    #=   println("\n")
    for i in N
        for j in N
            println("y[$i,$t] = ", round(sol.x[i,j], digits = 2), " ")
        end
        println("\n")
    end =#

    
    time_end = time_ns()
    
    stats.rf_time = (time_end - time_start) * 1e-9
    stats.rf_UB = sol.primal_bound

    # Check if solution is feasible (integer)
    integerSol = true
    N_INF= 0
    for j in N, i in N
        if value(sol.x[i,j]) > 0.01 && value(sol.x[i,j]) < 0.99
            integerSol = false
            N_INF +=1
            #println("FFAAIILL!!")
           # break
        end
    end
    #println("N_INF",N_INF)
    if integerSol == true && stats.sol_status != "Fail: !Int_TL" && stop_by_time_limit == false  
        stats.sol_status = "Success"
            println("\nRF found an integer solution with value: ", sol.primal_bound)
    else
        if stop_by_time_limit == false
            stats.sol_status = "Fail: !Int"
            println("\nRF could not find an integer solution")
        else
            stats.sol_status = "Fail: !Int_TL"
            println("\nRF stopped by time limit and could not find an integer solution")
        end
    end
    #= for i in N
        println("model[:u][i] = ",  value(model[:u][i]))
    end =#
    println("\n===================== Finished RF! ======================")

    
end




end #module
