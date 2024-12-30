module BestImprovement

using JuMP
using Data
using OutputStatistics



export bestImprovementSwap!, bestImprovementReinsertion!, bestImprovementOrOpt2!, bestImprovementOrOpt3!, bestImprovement2Opt!

function bestImprovementSwap!(data::InstanceData, data_type::InstanceTypeData, sol::StdFormModelSolution, stats::StatisticsData)

   # println(" Busca Local Swap ...")


    n = data_type.DIMENSION
    # é quantas cidades tem no total, o tamanho da rota é n+1
    
    rota_inicial = sol.route
    #println("rota_inicial: ", rota_inicial)

    

    bestDelta = 0
    besti = 0 # posicao inicial
    bestj = 0 # posicao final


    for idxi in 2:n - 2
        i = rota_inicial[idxi]

        for idxj in idxi+2:n

            j = rota_inicial[idxj]

            delta = 0
            
           # println("")
            delta -= data.c[rota_inicial[idxi-1],rota_inicial[idxi]]
            delta -= data.c[rota_inicial[idxi],rota_inicial[idxi+1]]
            delta -= data.c[rota_inicial[idxj-1],rota_inicial[idxj]]
            delta -= data.c[rota_inicial[idxj],rota_inicial[idxj+1]]
            delta += data.c[rota_inicial[idxi-1],rota_inicial[idxj]]
            delta += data.c[rota_inicial[idxj],rota_inicial[idxi+1]]
            delta += data.c[rota_inicial[idxj-1],rota_inicial[idxi]]
            delta += data.c[rota_inicial[idxi],rota_inicial[idxj+1]]
            
            #println(" delta[$i,$j] calculado: ", delta)
            
            if delta < bestDelta
                bestDelta = delta
                besti = i
                bestj = j
                #println("besti: ", besti)
                #println("bestj: ", bestj)
                #println("bestDelta: ", bestDelta)
            end
        end
       
        
    end
    for idxi in 2:n - 1
        i = rota_inicial[idxi]

        for idxj in idxi+1

            j = rota_inicial[idxj]

            delta = 0
            
           # println("")
            delta -= data.c[rota_inicial[idxi-1],rota_inicial[idxi]]
            delta -= data.c[rota_inicial[idxj],rota_inicial[idxj+1]]
            delta += data.c[rota_inicial[idxi-1],rota_inicial[idxj]]
            delta += data.c[rota_inicial[idxi],rota_inicial[idxj+1]]
            
            #println(" delta[$i,$j] calculado: ", delta)
            
            if delta < bestDelta
                bestDelta = delta
                besti = i
                bestj = j
               # println("besti: ", besti)
               # println("bestj: ", bestj)
                #println("bestDelta: ", bestDelta)
            end
        end
       
        
    end
   

    if bestDelta < 0

        # Remove a cidade na posição besti e insere na posição bestj
        pos1 = findfirst(==(besti), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        deleteat!(rota_inicial, pos1)  # Remove a cidade inicial do vetor cidade
        
        pos2 = findfirst(==(bestj), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        deleteat!(rota_inicial, pos2)  # Remove a cidade inicial do vetor cidade

        insert!(rota_inicial, pos2 , besti)  # Insere a cidade mais próxima após cidade_inicial

        insert!(rota_inicial, pos1 , bestj)  # Insere a cidade mais próxima após cidade_inicial
        #pos = findfirst(==(besti), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        #deleteat!(rota_inicial, pos-1)  # Remove a cidade inicial do vetor cidade
        
        
        
        
   

        sol.route = rota_inicial
        if stats.improvement == 0
            valor_inicial = stats.UB_HC
        else
            valor_inicial = stats.UB_LS
        end
        sol.primal_bound = sol.primal_bound + bestDelta 
        
        stats.UB_LS = sol.primal_bound
        stats.improvement = ((valor_inicial - stats.UB_LS)/valor_inicial)*100
       # println("besti: ", besti)
       # println("bestj: ", bestj)
       # println("best_UB_LS: ", stats.UB_LS)
       # println("route: ", sol.route)
        #println("")
    else
      #  println("Não houve melhora")
       # println("")
        stats.improvement = 0    
       # stats.UB_LS = sol.primal_bound
    end
    

    
    


end
function bestImprovement2Opt!(data::InstanceData, data_type::InstanceTypeData, sol::StdFormModelSolution, stats::StatisticsData)

   # println(" Busca Local heuristcs 2Opt ...")


    n = data_type.DIMENSION
    # é quantas cidades tem no total, o tamanho da rota é n+1
    
    rota_inicial = sol.route
    #println("rota_inicial: ", rota_inicial)

    

    bestDelta = 0
    besti = 0 # posicao inicial
    bestj = 0 # posicao final


    for idxi in 2:n -2
        i = rota_inicial[idxi]

        for idxj in idxi+2:n

            j = rota_inicial[idxj]

            delta = 0
            delta -= data.c[rota_inicial[idxi-1],rota_inicial[idxi]]
            delta -= data.c[rota_inicial[idxj],rota_inicial[idxj+1]]
            delta += data.c[rota_inicial[idxi-1],rota_inicial[idxj]]
            delta += data.c[rota_inicial[idxi],rota_inicial[idxj+1]]
            #println("delta: ", delta)
            if delta < bestDelta
                bestDelta = delta
                besti = i
                bestj = j
                #println("besti: ", besti)
                #println("bestj: ", bestj)
                #println("bestDelta: ", bestDelta)
            end
        end
       
        
    end
   
    if bestDelta < 0

        
        
        pos1 = findfirst(==(besti), rota_inicial)
        pos2 = findfirst(==(bestj), rota_inicial)
        a = reverse(rota_inicial[pos1:pos2])
        

        #println("a: ", a)
           #reverse!()
        for i=1:pos2-pos1

            deleteat!(rota_inicial, pos1)
        end
        #println("rota_inicial após remoção: ", rota_inicial)
        splice!(rota_inicial, pos1, a)
        
        

        sol.route = rota_inicial
        if stats.improvement == 0
            valor_inicial = stats.UB_HC
        else
            valor_inicial = stats.UB_LS
        end
        sol.primal_bound = sol.primal_bound + bestDelta 
        stats.UB_LS = sol.primal_bound
        stats.improvement = ((valor_inicial - stats.UB_LS)/valor_inicial)*100
        #println("best_UB_LS: ", stats.UB_LS)
       # println("besti: ", besti)
       # println("bestj: ", bestj)
        #println("route: ", sol.route)
       # println("")
       
    else
       # println("Não houve melhora")
       # println("")
        stats.improvement = 0    
    end
    

    
    


end
function bestImprovementReinsertion!(data::InstanceData, data_type::InstanceTypeData, sol::StdFormModelSolution, stats::StatisticsData)

   # println(" Busca Local Reinsertion ...")


    n = data_type.DIMENSION
    # é quantas cidades tem no total, o tamanho da rota é n+1
    
    rota_inicial = sol.route
    #println("rota_inicial: ", rota_inicial)

    

    bestDelta = 0
    besti = 0 # posicao inicial
    bestj = 0 # posicao final


    for idxi in 2:n - 2
        i = rota_inicial[idxi]

        for idxj in idxi+2:n

            j = rota_inicial[idxj]

            delta = 0
            delta -= data.c[rota_inicial[idxi-1],rota_inicial[idxi]]
            delta -= data.c[rota_inicial[idxi],rota_inicial[idxi+1]]
            delta -= data.c[rota_inicial[idxj],rota_inicial[idxj+1]]
            delta += data.c[rota_inicial[idxi-1],rota_inicial[idxi+1]]
            delta += data.c[rota_inicial[idxj],rota_inicial[idxi]]
            delta += data.c[rota_inicial[idxi],rota_inicial[idxj+1]]
            #println("delta: ", delta)
            if delta < bestDelta
                bestDelta = delta
                besti = i
                bestj = j
               #println("besti: ", besti)
               #println("bestj: ", bestj)
               #println("bestDelta: ", bestDelta)
            end
        end
       
        
    end
    for idxi in 2:n - 1
        i = rota_inicial[idxi]

        for idxj in idxi+1

            j = rota_inicial[idxj]

            delta = 0
            delta -= data.c[rota_inicial[idxi-1],rota_inicial[idxi]]
            delta -= data.c[rota_inicial[idxj],rota_inicial[idxj+1]]
            delta += data.c[rota_inicial[idxi-1],rota_inicial[idxj]]
            delta += data.c[rota_inicial[idxi],rota_inicial[idxj+1]]
            #println("delta: ", delta)
            if delta < bestDelta
                bestDelta = delta
                besti = i
                bestj = j
               #println("besti: ", besti)
               #println("bestj: ", bestj)
               #println("bestDelta: ", bestDelta)
            end
        end
       
        
    end
   
    if bestDelta < 0

        # Remove a cidade na posição besti e insere na posição bestj
        pos1 = findfirst(==(besti), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        deleteat!(rota_inicial, pos1)  # Remove a cidade inicial do vetor cidade

        pos2 = findfirst(==(bestj), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        
        insert!(rota_inicial, pos2+1 , besti)  # Insere a cidade mais próxima após cidade_inicial

        sol.route = rota_inicial
        if stats.improvement == 0
            valor_inicial = stats.UB_HC
        else
            valor_inicial = stats.UB_LS
        end
        sol.primal_bound = sol.primal_bound + bestDelta 
        stats.UB_LS = sol.primal_bound
        stats.improvement = ((valor_inicial - stats.UB_LS)/valor_inicial)*100
       # println("best_UB_LS: ", stats.UB_LS)
        #println("besti: ", besti)
       # println("bestj: ", bestj)
        #println("route: ", sol.route)
       # println("")
    else
       # println("Não houve melhora")
       # println("")
        stats.improvement = 0   
       # stats.UB_LS = sol.primal_bound 
    end
    

    
    


end
function bestImprovementOrOpt2!(data::InstanceData, data_type::InstanceTypeData, sol::StdFormModelSolution, stats::StatisticsData)

   # println(" Busca Local OrOpt2 ...")


    n = data_type.DIMENSION
    # é quantas cidades tem no total, o tamanho da rota é n+1
    
    rota_inicial = sol.route
    #println("rota_inicial: ", rota_inicial)

    

    bestDelta = 0
    besti1 = 0 # posicao inicial
    besti2 = 0 # posicao inicial
    bestj = 0 # posicao final


    for idxi in 2:n - 2
        i1 = rota_inicial[idxi]
        i2 = rota_inicial[idxi+1]
        for idxj in idxi+2:n

            j = rota_inicial[idxj]

            delta = 0
            delta -= data.c[rota_inicial[idxi-1],rota_inicial[idxi]]
            delta -= data.c[rota_inicial[idxi+1],rota_inicial[idxi+2]]
            delta -= data.c[rota_inicial[idxj],rota_inicial[idxj+1]]
            delta += data.c[rota_inicial[idxi-1],rota_inicial[idxi+2]]
            delta += data.c[rota_inicial[idxj],rota_inicial[idxi]]
            delta += data.c[rota_inicial[idxi+1],rota_inicial[idxj+1]]
            #println("delta: ", delta)
            if delta < bestDelta
                bestDelta = delta
                besti1 = i1
                besti2 = i2
                bestj = j
               # println("besti: ", besti1)
               # println("besti: ", besti2)
               # println("bestj: ", bestj)
               # println("bestDelta: ", bestDelta)
            end
        end
       
        
    end
    
    if bestDelta < 0

        # Remove a cidade na posição besti e insere na posição bestj
        pos = findfirst(==(besti1), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        deleteat!(rota_inicial, pos)  # Remove a cidade inicial do vetor cidade
        pos = findfirst(==(besti2), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        deleteat!(rota_inicial, pos)  # Remove a cidade inicial do vetor cidade

        pos = findfirst(==(bestj), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        insert!(rota_inicial, pos + 1, besti1)  # Insere a cidade mais próxima após cidade_inicial
        pos = findfirst(==(bestj), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        insert!(rota_inicial, pos + 2, besti2)  # Insere a cidade mais próxima após cidade_inicial



        sol.route = rota_inicial

        if stats.improvement == 0
            valor_inicial = stats.UB_HC
        else
            valor_inicial = stats.UB_LS
        end
        sol.primal_bound = sol.primal_bound + bestDelta 
        stats.UB_LS = sol.primal_bound
        stats.improvement = ((valor_inicial - stats.UB_LS)/valor_inicial)*100
       # println("besti: ", besti1)
       # println("bestj: ", bestj)
       # println("best_UB_LS: ", stats.UB_LS)
       # println("route: ", sol.route)
      #  println("")
       
    else
       # println("Não houve melhora")
      #  println("")
        stats.improvement = 0    
        
    end
    

    
    


end
function bestImprovementOrOpt3!(data::InstanceData, data_type::InstanceTypeData, sol::StdFormModelSolution, stats::StatisticsData)

   # println(" Busca Local  OrOpt3 ...")


    n = data_type.DIMENSION
    # é quantas cidades tem no total, o tamanho da rota é n+1
    
    rota_inicial = sol.route
   # println("rota_inicial: ", rota_inicial)

    

    bestDelta = 0
    besti1 = 0 # posicao inicial
    besti2 = 0 # posicao inicial
    besti3 = 0 # posicao inicial
    bestj = 0 # posicao final


    for idxi in 2:n - 3
        i1 = rota_inicial[idxi]
        i2 = rota_inicial[idxi+1]
        i3 = rota_inicial[idxi+2]
        for idxj in idxi+3:n

            j = rota_inicial[idxj]

            delta = 0
            delta -= data.c[rota_inicial[idxi-1],rota_inicial[idxi]]
            delta -= data.c[rota_inicial[idxi+2],rota_inicial[idxi+3]]
            delta -= data.c[rota_inicial[idxj],rota_inicial[idxj+1]]
            delta += data.c[rota_inicial[idxi-1],rota_inicial[idxi+3]]
            delta += data.c[rota_inicial[idxj],rota_inicial[idxi]]
            delta += data.c[rota_inicial[idxi+2],rota_inicial[idxj+1]]
            #println("delta: ", delta)
            if delta < bestDelta
                bestDelta = delta
                besti1 = i1
                besti2 = i2
                besti3 = i3
                bestj = j
                #println("besti: ", besti1)
                #println("besti: ", besti2)
                #println("besti: ", besti3)
                #println("bestj: ", bestj)
                #println("bestDelta: ", bestDelta)
            end
        end
       
        
    end
   
    if bestDelta < 0

        # Remove a cidade na posição besti e insere na posição bestj
        pos = findfirst(==(besti1), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        deleteat!(rota_inicial, pos)  # Remove a cidade inicial do vetor cidade
        pos = findfirst(==(besti2), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        deleteat!(rota_inicial, pos)  # Remove a cidade inicial do vetor cidade
        pos = findfirst(==(besti3), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        deleteat!(rota_inicial, pos)  # Remove a cidade inicial do vetor cidade

        pos = findfirst(==(bestj), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        insert!(rota_inicial, pos + 1, besti1)  # Insere a cidade mais próxima após cidade_inicial
        pos = findfirst(==(bestj), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        insert!(rota_inicial, pos + 2, besti2)  # Insere a cidade mais próxima após cidade_inicial
        pos = findfirst(==(bestj), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
        insert!(rota_inicial, pos + 3, besti3)  # Insere a cidade mais próxima após cidade_inicial

        

        sol.route = rota_inicial

        if stats.improvement == 0
            valor_inicial = stats.UB_HC
        else
            valor_inicial = stats.UB_LS
        end
        sol.primal_bound = sol.primal_bound + bestDelta 
        stats.UB_LS = sol.primal_bound
      
        stats.improvement = ((valor_inicial - stats.UB_LS)/valor_inicial)*100
       # println("best_UB_LS: ", stats.UB_LS)
       # println("besti: ", besti1)
        #println("bestj: ", bestj)
        #println("route: ", sol.route)
       # println("")
       
    else
       # println("Não houve melhora")
       # println("")
        stats.improvement = 0
        
    end
    

    
    


end

end # module