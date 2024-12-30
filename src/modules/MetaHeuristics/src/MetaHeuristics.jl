module MetaHeuristics

using JuMP
using Data
using OutputStatistics
using BestImprovement
using Random
using Disturbance
using Constructive



export ILS!, MultiStart!, SA!, VNS!, GRASP!, TABU!

function ILS!(data::InstanceData, data_type::InstanceTypeData, sol::StdFormModelSolution, stats::StatisticsData)

    Constructive.solve_IMB!(data, data_type, sol, stats)

    

    n = data_type.DIMENSION
    # é quantas cidades tem no total, o tamanho da rota é n+1
    
    #rota_atua = sol.route

    best_rota = sol.route
    best_UB = stats.UB_HC
    iterIls = 0
    if n >= 150
        maxIterIls = div(n, 2)
    else
        maxIterIls = n
    end
   
    while iterIls <= maxIterIls


        maxIter = 50
        Iter = 0

        while Iter <= maxIter
            
            
            RD = rand(1:5)
            if RD == 1
                BestImprovement.bestImprovementSwap!(data, data_type, sol, stats)
            elseif RD == 2
                BestImprovement.bestImprovementReinsertion!(data, data_type, sol, stats)
            elseif RD == 3
                BestImprovement.bestImprovementOrOpt2!(data, data_type, sol, stats)
            elseif RD == 4
            BestImprovement.bestImprovementOrOpt3!(data, data_type, sol, stats)
            elseif RD == 5
                BestImprovement.bestImprovement2Opt!(data, data_type, sol, stats)
            end

            if stats.improvement == 0
                Iter += 1
            else
                Iter = 0
            end
            #println("stats.UB_LS aqui !!!!", stats.UB_LS)

          
            

        end
        
        
        solution_route = sol.route
        #println("best_UB = ", best_UB) 

        if stats.UB_LS < best_UB
            best_UB = stats.UB_LS
            best_rota = solution_route
            iterIls = 0


        else
            stats.UB_LS = best_UB
            sol.route = best_rota
            iterIls += 1
        end
        #println("new best_UB = ", best_UB)
        if iterIls <= maxIterIls
            Disturbance.DoubleBridge!(data, data_type, sol, stats)
        end
        

        
       
        
        
    end

   
        
    stats.UB_LS = sol.primal_bound
    println("solução !!!!", stats.UB_LS)
    stats.improvement = (( stats.UB_HC - stats.UB_LS)/ stats.UB_HC)*100

end
function MultiStart!(data::InstanceData, data_type::InstanceTypeData, sol::StdFormModelSolution, stats::StatisticsData)

    Constructive.solve_VMP!(data, data_type, sol, stats)

    

    n = data_type.DIMENSION
    # é quantas cidades tem no total, o tamanho da rota é n+1
    
    #rota_atua = sol.route

    best_rota = sol.route
    best_UB = stats.UB_HC
    iterMultiStart = 0
    if n >= 150
        maxIterMultiStart = div(n, 2)
    else
        maxIterMultiStart = n
    end
     
    while iterMultiStart <= maxIterMultiStart

        Constructive.solve_random!(data, data_type, sol, stats)

        





        maxIter = 10
        Iter = 0

        while Iter <= maxIter
            
            
            RD = rand(1:5)
            if RD == 1
                BestImprovement.bestImprovementSwap!(data, data_type, sol, stats)
            elseif RD == 2
                BestImprovement.bestImprovementReinsertion!(data, data_type, sol, stats)
            elseif RD == 3
                BestImprovement.bestImprovementOrOpt2!(data, data_type, sol, stats)
            elseif RD == 4
                BestImprovement.bestImprovementOrOpt3!(data, data_type, sol, stats)
            elseif RD == 5
                BestImprovement.bestImprovement2Opt!(data, data_type, sol, stats)
            end

            if stats.improvement == 0
                Iter += 1
            else
                Iter = 0
            end
            #println("stats.UB_LS aqui !!!!", stats.UB_LS)

            if stats.UB_LS !=  -1e8
                if stats.UB_LS < 0 #12
                   # println("parar")
                    break
                end 
            end
            

        end
       
        solution_route = sol.route
       # println("best_UB = ", best_UB) 

        if stats.UB_LS < best_UB
            best_UB = stats.UB_LS
            best_rota = solution_route
            iterMultiStart = 0


        else
            #= stats.UB_LS = best_UB
            sol.route = best_rota =#
            iterMultiStart += 1
        end
      #  println("new best_UB = ", best_UB)
     
        

        
       
        
        
    end

   
        
    stats.UB_LS = sol.primal_bound
    stats.improvement = (( stats.UB_HC - stats.UB_LS)/ stats.UB_HC)*100
    println("solução !!!!", stats.UB_LS)
end
function SA!(data::InstanceData, data_type::InstanceTypeData, sol::StdFormModelSolution, stats::StatisticsData)

    #println(" contruir solução incial SA ")
    RD = 1 #rand(1:4)

        if RD == 1
            Constructive.solve_VMP!(data, data_type, sol, stats)
        elseif RD == 2
            Constructive.solve_BN!(data, data_type, sol, stats)
        elseif RD == 3
            Constructive.solve_IMB!(data, data_type, sol, stats)
        elseif RD == 4
            Constructive.solve_random!(data, data_type, sol, stats)
        end 

    n = data_type.DIMENSION
    # é quantas cidades tem no total, o tamanho da rota é n+1
    
    #rota_atua = sol.route

    best_rota = sol.route
    best_UB = stats.UB_HC
    iterSA = 0
    if n >= 150
        maxIterSA = div(n, 2)
    else
        maxIterSA = n
    end
    
        #println("busca local SA")
    T = 10000
    while T > 0.01
        while iterSA <= maxIterSA

            iterSA += 1
            RD = rand(1:5)

            if RD == 1
                BestImprovement.bestImprovementSwap!(data, data_type, sol, stats)
            elseif RD == 2
                BestImprovement.bestImprovementReinsertion!(data, data_type, sol, stats)
            elseif RD == 3
                BestImprovement.bestImprovementOrOpt2!(data, data_type, sol, stats)
            elseif RD == 4
                BestImprovement.bestImprovementOrOpt3!(data, data_type, sol, stats)
            elseif RD == 5
                BestImprovement.bestImprovement2Opt!(data, data_type, sol, stats)
            end
           
            delta = best_UB - stats.UB_LS

            if delta > 0
                
                best_UB = stats.UB_LS
                best_rota = sol.route
                

            else
                if rand() < exp(delta/T)
                    best_UB = stats.UB_LS
                    best_rota = sol.route
                end
            end
            
        end

        T = 0.95    *T

    end
        
    stats.UB_LS = sol.primal_bound
    stats.improvement = (( stats.UB_HC - stats.UB_LS)/ stats.UB_HC)*100

end
function VNS!(data::InstanceData, data_type::InstanceTypeData, sol::StdFormModelSolution, stats::StatisticsData)

    Constructive.solve_IMB!(data, data_type, sol, stats)

    

    n = data_type.DIMENSION
    # é quantas cidades tem no total, o tamanho da rota é n+1
    
    #rota_atua = sol.route

    best_rota = sol.route
    best_UB = stats.UB_HC
    iterVNS = 0
    
    maxIterVNS = 50
    while iterVNS <= maxIterVNS


        maxIter = 10
        Iter = 0

        RD = 1
        while Iter <= maxIter
            
            
            RD 
            if RD == 1
                BestImprovement.bestImprovementSwap!(data, data_type, sol, stats)
            elseif RD == 2
                BestImprovement.bestImprovementReinsertion!(data, data_type, sol, stats)
            elseif RD == 3
                BestImprovement.bestImprovementOrOpt2!(data, data_type, sol, stats)
            elseif RD == 4
                BestImprovement.bestImprovementOrOpt3!(data, data_type, sol, stats)
            elseif RD == 5
                BestImprovement.bestImprovement2Opt!(data, data_type, sol, stats)
            end
                RD += 1
                if RD > 5
                    RD = 1
                end
            if stats.improvement == 0
                Iter += 1
            else
                Iter = 0
            end
            #println("stats.UB_LS aqui !!!!", stats.UB_LS)

            if stats.UB_LS !=  -1e8
                if stats.UB_LS < 0 #12
                   # println("parar")
                    break
                end 
            end
            

        end
        if stats.UB_LS !=  -1e8
            if stats.UB_LS < 0 #12
               # println("parar")
                break
            end 
        end
        
        solution_route = sol.route
        println("best_UB = ", best_UB) 

        if stats.UB_LS < best_UB
            best_UB = stats.UB_LS
            best_rota = solution_route
            iterVNS = 0


        else
            stats.UB_LS = best_UB
            sol.route = best_rota
            iterVNS += 1
        end
        println("new best_UB = ", best_UB)
       
        

        
       
        
        
    end

   
        
    stats.UB_LS = sol.primal_bound
    stats.improvement = (( stats.UB_HC - stats.UB_LS)/ stats.UB_HC)*100

end
function GRASP!(data::InstanceData, data_type::InstanceTypeData, sol::StdFormModelSolution, stats::StatisticsData)

    
 
     
 
    n = data_type.DIMENSION
    # é quantas cidades tem no total, o tamanho da rota é n+1
     
    #rota_atua = sol.route
 
    # best_rota = sol.route
    #best_UB = stats.UB_HC
    iterGRASP = 0
    best_UB = Inf
    best_rota = []
    if n >= 150
        maxIterGRASP = div(n, 2)
    else
        maxIterGRASP = n
    end
    
    while iterGRASP <= maxIterGRASP
 
        Constructive.solve_C_GRASP!(data, data_type, sol, stats)
        #println("solução inicial GRASP",stats.UB_HC)
        f = stats.UB_HC
        s = []
        improvement = true
        while improvement == true 
            f1 = Inf
            f2 = Inf
            f3 = Inf
            f4 = Inf
            f5 = Inf
            s1 = []
            s2 = []
            s3 = []
            s4 = []
            s5 = []
            

            BestImprovement.bestImprovementSwap!(data, data_type, sol, stats)
            if stats.UB_LS < f
                f1 = stats.UB_LS
                s1 = sol.route
               # println("f1  !!!!", stats.UB_LS)
            end
            
            BestImprovement.bestImprovementReinsertion!(data, data_type, sol, stats)
            if stats.UB_LS < f
                f2 = stats.UB_LS
                s2 = sol.route
              #  println("f2  !!!!", stats.UB_LS)
            end
            
            BestImprovement.bestImprovementOrOpt2!(data, data_type, sol, stats)
            if stats.UB_LS < f
                f3 = stats.UB_LS
                s3 = sol.route
               # println("f3  !!!!", stats.UB_LS)
            end

            BestImprovement.bestImprovementOrOpt3!(data, data_type, sol, stats)
            if stats.UB_LS < f
                f4 = stats.UB_LS
                s4 = sol.route
                #println("f4  !!!!", stats.UB_LS)
            end
          
            BestImprovement.bestImprovement2Opt!(data, data_type, sol, stats)
            if stats.UB_LS < f
                f5 = stats.UB_LS
                s5 = sol.route
               # println("f5  !!!!", stats.UB_LS)
            end
           
            f_values = [f1, f2, f3, f4, f5]
            s_values = [s1, s2, s3, s4, s5]

            # Filter out Inf values
            filtered_f_values = filter(x -> x != Inf, f_values)
            filtered_s_values = [s_values[i] for i in 1:length(f_values) if f_values[i] != Inf]

           # println("filtered_f_values  !!!!", filtered_f_values)

            if !isempty(filtered_f_values)
                f, index = findmin(filtered_f_values)
               # println("f  !!!!", f)
                s = filtered_s_values[index]
               # println("s  !!!!", s)
                improvement = true
            else
                improvement = false
            end

            
            
            #println("stats.UB_LS aqui !!!!", stats.UB_LS)
            
            
            
            
            
           
             
        end
         
        iterGRASP += 1    
         
        solution_route = sol.route
       # println("best_UB = ", best_UB) 
 
        if f < best_UB
            best_UB = f
            best_rota = s
           
        else
            stats.UB_LS = best_UB
            sol.route = best_rota
            iterGRASP += 1
        end
       # println("new best_UB = ", best_UB)
         
                  
         
    end
 
    
         
    stats.UB_LS = sol.primal_bound
    stats.improvement = (( stats.UB_HC - stats.UB_LS)/ stats.UB_HC)*100
 
 
end 
function TABU!(data::InstanceData, data_type::InstanceTypeData, sol::StdFormModelSolution, stats::StatisticsData)

    Constructive.solve_IMB!(data, data_type, sol, stats)
 
     
 
    n = data_type.DIMENSION
    # é quantas cidades tem no total, o tamanho da rota é n+1
     
 
    best_rota = sol.route
    best_UB = stats.UB_HC
    iterTABU = 0
    maxIterTABU = 50
    bestTABU = 0
    TABU_Liste = []
    T_TABU_Liste = 2

    while iterTABU - bestTABU <= maxIterTABU 

        f1 = Inf
        f2 = Inf
        f3 = Inf
        f4 = Inf
        f5 = Inf
        s1 = []
        s2 = []
        s3 = []
        s4 = []
        s5 = []
        
        
        

            
        if !(1 in TABU_Liste)
            BestImprovement.bestImprovementSwap!(data, data_type, sol, stats)
            f1 = stats.UB_LS
            s1 = sol.route
            # println("f1  !!!!", stats.UB_LS)
        end
        
        
        if !(2 in TABU_Liste)
            BestImprovement.bestImprovementReinsertion!(data, data_type, sol, stats)
            f2 = stats.UB_LS
            s2 = sol.route
            #  println("f2  !!!!", stats.UB_LS)
        end
        
       
        if !(3 in TABU_Liste)
            BestImprovement.bestImprovementOrOpt2!(data, data_type, sol, stats)
            f3 = stats.UB_LS
            s3 = sol.route
            # println("f3  !!!!", stats.UB_LS)
        end

       
        if !(4 in TABU_Liste)
            BestImprovement.bestImprovementOrOpt3!(data, data_type, sol, stats)
            f4 = stats.UB_LS
            s4 = sol.route
            #println("f4  !!!!", stats.UB_LS)
        end
        
       
        if !(5 in TABU_Liste)
            BestImprovement.bestImprovement2Opt!(data, data_type, sol, stats)
            f5 = stats.UB_LS
            s5 = sol.route
            # println("f5  !!!!", stats.UB_LS)
        end
        
        f_values = [f1, f2, f3, f4, f5]
        s_values = [s1, s2, s3, s4, s5]

        # Filter out Inf values
        filtered_f_values = filter(x -> x != Inf, f_values)
        filtered_s_values = [s_values[i] for i in 1:length(f_values) if f_values[i] != Inf]

        # println("filtered_f_values  !!!!", filtered_f_values

        
        f, index = findmin(filtered_f_values)
        # println("f  !!!!", f)
        s = filtered_s_values[index]

        push!(TABU_Liste, index)
        if length(TABU_Liste) > T_TABU_Liste
            popfirst!(TABU_Liste)
        end
        
             
         
       # solution_route = sol.route
       # println("best_UB = ", best_UB) 
 
        if f < best_UB
            best_UB = f
            best_rota = s
            bestTABU = iterTABU
            #= else
            stats.UB_LS = best_UB
            sol.route = best_rota
             =#
        end
       # println("new best_UB = ", best_UB)
         
                
        iterTABU += 1    
       
    end
 
    
         
    stats.UB_LS = sol.primal_bound
    stats.improvement = (( stats.UB_HC - stats.UB_LS)/ stats.UB_HC)*100
 
    
end
 
end# module