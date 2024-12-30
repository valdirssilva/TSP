module Constructive

using JuMP
using Data
using OutputStatistics
using Random


export solve_VMP!, solve_BN!, solve_IMB!, solver_random!, solve_C_GRASP!



function solve_VMP!(data::InstanceData,  data_type::InstanceTypeData,sol::StdFormModelSolution, stats::StatisticsData)

  #  println(" Construtive heuristcs nearest neighbor...")


    n = data_type.DIMENSION
    cidades_visitar = collect(1:n)  # Cria um vetor com N itens que vão de 1 até N
    #sol.x = zeros(Int, n, n) # Cria uma matriz de zeros com N linhas e N colunas
    
    cidade_inicial = first(cidades_visitar)  # Pega o primeiro valor do vetor cidade
    cidades_visitadas = [cidade_inicial]
    

    distancia_total = 0
    cidade_inicial = 1
    while length(cidades_visitar) > 0

        cidade_mais_proxima = 0
        
        i = cidade_inicial
        #println("i",i)
        distancia_mais_proxima = Inf

        if length(cidades_visitar) != 1

            for j in cidades_visitar #first(cidade):last(cidade)
                if j != i
                    distancia = data.c[i, j]
                    #println("data.c[$i, $j]",data.c[i, j])
                    if distancia < distancia_mais_proxima
                        distancia_mais_proxima = distancia
                        cidade_mais_proxima = j
                        
                    end
                end
            end
        else

            distancia_mais_proxima = data.c[i, 1]
            cidade_mais_proxima = 1
        end




        distancia_total = distancia_total + distancia_mais_proxima
       # println("Distancia total: ", distancia_total)

       if length(cidades_visitadas) > 0

            pos = findfirst(==(cidade_inicial), cidades_visitar)  # Encontra a posição do valor cidade_inicial no vetor cidade
            # println("Posição: ", pos)

            deleteat!(cidades_visitar, pos)  # Remove a cidade inicial do vetor cidade
        # println("Cidades: ", cidade)

       end
       if length(cidades_visitadas) > 0
            pos = findfirst(==(cidade_inicial), cidades_visitadas)  # Encontra a posição do valor cidade_inicial no vetor cidade
            # println("Posição: ", pos)
            if pos > 1
                insert!(cidades_visitadas, pos + 1, cidade_mais_proxima)  # Insere a cidade mais próxima após cidade_inicial
            else
                insert!(cidades_visitadas, pos, cidade_mais_proxima)  # Insere a cidade mais próxima após cidade_inicial
            end
        end
    
        sol.x[cidade_inicial, cidade_mais_proxima] = 1
        cidade_inicial = cidade_mais_proxima

        #println("cidades_visitar", cidades_visitar)
        #println("cidades_visitadas", cidades_visitadas)


    end
    sol.route = cidades_visitadas
    sol.primal_bound = distancia_total
    stats.UB_HC = sol.primal_bound
   # println("best_UB: ", stats.UB_HC)
    return
end

function solve_BN!(data::InstanceData,  data_type::InstanceTypeData,sol::StdFormModelSolution, stats::StatisticsData)

    #println(" Construtive heuristcs Bellmore e Nemhauser...")


    n = data_type.DIMENSION
    cidades_visitar = collect(1:n)  # Cria um vetor com N itens que vão de 1 até N
    
    sol.x = zeros(Int, n, n) # Cria uma matriz de zeros com N linhas e N colunas
    
    cidade_inicial = first(cidades_visitar)  # Pega o primeiro valor do vetor cidade
    cidades_visitadas = [cidade_inicial]
    deleteat!(cidades_visitar, 1)


    while length(cidades_visitar) > 0
        cidade_mais_proxima = 0
        
        
        #println("i",i)
        distancia_mais_proxima = Inf 

       # if length(cidades_visitar) != 0
            for i in [first(cidades_visitadas), last(cidades_visitadas)]
                for j in cidades_visitar #first(cidade):last(cidade)
                    if j != i
                        distancia = data.c[i, j]
                      # println("data.c[$i, $j]",data.c[i, j])
                        if distancia < distancia_mais_proxima
                            distancia_mais_proxima = distancia
                            cidade_mais_proxima = j
                            cidade_inicial = i
                           #println("Cidade inicial: ", cidade_inicial, " Cidade mais proxima: ", cidade_mais_proxima)
                            
                        end
                    end
                end
            end

        #end
        #println("distancia mais proxima: ", distancia_mais_proxima)
        pos = findfirst(==(cidade_inicial), cidades_visitadas)  # Encontra a posição do valor cidade_inicial no vetor cidade
           
        if pos > 1
            insert!(cidades_visitadas, pos +1, cidade_mais_proxima)  # Insere a cidade mais próxima após cidade_inicial
        else
            insert!(cidades_visitadas, pos , cidade_mais_proxima)  # Insere a cidade mais próxima após cidade_inicial
        end
    

       if length(cidades_visitar) > 0
            pos = findfirst(==(cidade_mais_proxima), cidades_visitar)  # Encontra a posição do valor cidade_inicial no vetor cidade
           # println("Posição: ", pos)
        
            deleteat!(cidades_visitar, pos)  # Remove a cidade inicial do vetor cidade
       end
       

       



       #println("Cidade inicial: ", cidade_inicial, " Cidade mais proxima: ", cidade_mais_proxima)
       #println("Cidades visitadas: ", cidades_visitadas)
       #println("Cidades a visitar: ", cidades_visitar)

       cidade_inicial = cidade_mais_proxima
       # println("Cidades: ", cidade)
        sol.x[cidade_inicial, cidade_mais_proxima] = 1
        


    end

    
    
        fechar_ciclo = first(cidades_visitadas)
        insert!(cidades_visitadas, length(cidades_visitadas) + 1, fechar_ciclo)  # Insere a cidade mais próxima após cidade_inicial

    
    #cidades_visitadas [1 2 4 5 6 8 7  2 6 4 ]
    distancia_total = 0
    for i in 1:(length(cidades_visitadas) - 1)
        distancia_total += data.c[cidades_visitadas[i], cidades_visitadas[i + 1]]
    end
    sol.x = zeros(Int, n, n) # Cria uma matriz de zeros com N linhas e N colunas
    for i in 1:(length(cidades_visitadas) - 1)
        sol.x[cidades_visitadas[i], cidades_visitadas[i + 1]] = 1
       # println("cidade inicial: ", cidades_visitadas[i], " cidade visitada: ", cidades_visitadas[i + 1]) 
       # println("sol[",cidades_visitadas[i],",",cidades_visitadas[i + 1],"] = ",sol.x[cidades_visitadas[i], cidades_visitadas[i + 1]])
       
    end
    
    
   # println("Cidades visitadas: ", cidades_visitadas)
    sol.route = cidades_visitadas
    sol.primal_bound = distancia_total
    stats.UB_HC = sol.primal_bound
   # println("best_UB: ", stats.UB_HC)
    return
end
function solve_IMB!(data::InstanceData,  data_type::InstanceTypeData,sol::StdFormModelSolution, stats::StatisticsData)

   # println(" Construtive heuristcs Inserção Mais Barata...")

   # println(" etapa 1")
    n = data_type.DIMENSION
    cidades_visitar = collect(1:n)  # Cria um vetor com N itens que vão de 1 até N
    
    sol.x = zeros(Int, n, n) # Cria uma matriz de zeros com N linhas e N colunas
    
    cidade_inicial = first(cidades_visitar)  # Pega o primeiro valor do vetor cidade
    cidades_visitadas = [cidade_inicial]
    deleteat!(cidades_visitar, 1)


   
    cidade_mais_proxima = 0
        
        
    #println("i",i)
    distancia_mais_proxima = Inf

    if length(cidades_visitar) != 0
        for i in [first(cidades_visitadas), last(cidades_visitadas)]
            for j in cidades_visitar #first(cidade):last(cidade)
                if j != i
                    distancia = data.c[i, j]
                    # println("data.c[$i, $j]",data.c[i, j])
                    if distancia < distancia_mais_proxima
                        distancia_mais_proxima = distancia
                        cidade_mais_proxima = j
                        cidade_inicial = i
                        #println("Cidade inicial: ", cidade_inicial, " Cidade mais proxima: ", cidade_mais_proxima)
                        
                    end
                end
            end
        end

    end
    #println("distancia mais proxima: ", distancia_mais_proxima)

    if length(cidades_visitadas) > 0

        insert!(cidades_visitadas, 2, cidade_mais_proxima)  # Insere a cidade mais próxima após cidade_inicial
  
    end
    

    if length(cidades_visitar) > 0
        pos = findfirst(==(cidade_mais_proxima), cidades_visitar)  # Encontra a posição do valor cidade_inicial no vetor cidade
        # println("Posição: ", pos)
    
        deleteat!(cidades_visitar, pos)  # Remove a cidade inicial do vetor cidade
    end
   # println("cidades_visitar", cidades_visitar)
    #println("cidades_visitadas", cidades_visitadas)

   # println("estapa 2")
   
    distancia_mais_proxima = Inf
    

        for j in cidades_visitar
            
            distancia = data.c[first(cidades_visitadas), j] + data.c[j,last(cidades_visitadas)]
            #println("distancia de $(first(cidades_visitadas)), $j e $(last(cidades_visitadas)) = ", data.c[first(cidades_visitadas), j] + data.c[j,last(cidades_visitadas)])  
            
            if distancia < distancia_mais_proxima
                distancia_mais_proxima = distancia
                cidade_mais_proxima = j
                # println( " Cidade mais proxima: ", cidades_visitar[cidade_mais_proxima])
                
            end
            
        end


       # println(" Cidade mais proxima: ", cidade_mais_proxima)

        if length(cidades_visitadas) > 0

            insert!(cidades_visitadas,  3, cidade_mais_proxima)  # Insere a cidade mais próxima após cidade_inicial
        end
        
        if length(cidades_visitar) > 0
            pos = findfirst(==(cidade_mais_proxima), cidades_visitar)  # Encontra a posição do valor cidade_inicial no vetor cidade
            # println("Posição: ", pos)
        
            deleteat!(cidades_visitar, pos)  # Remove a cidade inicial do vetor cidade
        end
      #  println("cidades_visitar", cidades_visitar)
      #  println("cidades_visitadas", cidades_visitadas)

    
   
   # println("etapa 3")
    while length(cidades_visitar) > 0
    
        
        distancia_mais_proxima = Inf

        for idx in 1:(length(cidades_visitadas) - 1)
            i = cidades_visitadas[idx]
            j = cidades_visitadas[idx + 1]

            for k in cidades_visitar
                distancia = data.c[i, k] + data.c[ k,j] - data.c[i, j]
               # println("distancia de $i, $k e $j = ", distancia) 
                if distancia < distancia_mais_proxima
                    distancia_mais_proxima = distancia
                    cidade_mais_proxima = k
                    cidade_inicial = i
                    #println("Cidade inicial: ", cidade_inicial, " Cidade mais proxima: ", cidade_mais_proxima)
                end
            end
        end
        i = last(cidades_visitadas)
        j = first(cidades_visitadas)

        for k in cidades_visitar
            distancia = data.c[i, k] + data.c[ k,j] - data.c[i, j]
           # println("distancia de $i, $k e $j = ", distancia)
            if distancia < distancia_mais_proxima
                distancia_mais_proxima = distancia
                cidade_mais_proxima = k
                cidade_inicial = i
                #println("Cidade inicial: ", cidade_inicial, " Cidade mais proxima: ", cidade_mais_proxima)
            end
        end



       # println("Cidade inicial: ", cidade_inicial, " Cidade mais proxima: ", cidade_mais_proxima)
                        


        if length(cidades_visitadas) > 0
            pos = findfirst(==(cidade_inicial), cidades_visitadas)  # Encontra a posição do valor cidade_inicial no vetor cidade
            # println("Posição: ", pos)
            insert!(cidades_visitadas, pos + 1, cidade_mais_proxima)  # Insere a cidade mais próxima após cidade_inicial
       #=      if pos > 1
                insert!(cidades_visitadas, pos + 1, cidade_mais_proxima)  # Insere a cidade mais próxima após cidade_inicial
            else
                insert!(cidades_visitadas, pos, cidade_mais_proxima)  # Insere a cidade mais próxima após cidade_inicial
            end =#
            
       end
        
      # println("cidades_visitadas", cidades_visitadas)

       if length(cidades_visitar) > 0
            pos = findfirst(==(cidade_mais_proxima), cidades_visitar)  # Encontra a posição do valor cidade_inicial no vetor cidade
           # println("Posição: ", pos)
        
            deleteat!(cidades_visitar, pos)  # Remove a cidade inicial do vetor cidade
       end

       #println("cidades_visitar", cidades_visitar)
 




    end
   # println("etapa 4")

    
    
        
    insert!(cidades_visitadas, length(cidades_visitadas) + 1, 1)  # Insere a cidade mais próxima após cidade_inicial

    
    #cidades_visitadas [1 2 4 5 6 8 7  2 6 4 ]
    distancia_total = 0
    for i in 1:(length(cidades_visitadas) - 1)
        distancia_total += data.c[cidades_visitadas[i], cidades_visitadas[i + 1]]
    end 
    sol.x = zeros(Int, n, n) # Cria uma matriz de zeros com N linhas e N colunas
    for i in 1:(length(cidades_visitadas) - 1)
        sol.x[cidades_visitadas[i], cidades_visitadas[i + 1]] = 1
       # println("cidade inicial: ", cidades_visitadas[i], " cidade visitada: ", cidades_visitadas[i + 1]) 
        #println("sol[",cidades_visitadas[i],",",cidades_visitadas[i + 1],"] = ",sol.x[cidades_visitadas[i], cidades_visitadas[i + 1]])
       
    end
    sol.route = cidades_visitadas

      #  println("Cidades visitadas: ", cidades_visitadas)

    sol.primal_bound = distancia_total
    stats.UB_HC = sol.primal_bound
   # println("best_UB: ", stats.UB_HC)
   

    return
end
function solve_random!(data::InstanceData,  data_type::InstanceTypeData,sol::StdFormModelSolution, stats::StatisticsData)

    #println(" Construtive heuristcs random...")


    n = data_type.DIMENSION
    cidades_visitar = collect(1:n)  # Cria um vetor com N itens que vão de 1 até N
    cidades_visitadas = shuffle(cidades_visitar)
    
    #println("cidades visitar: ", cidades_visitar)
   # println("cidades visitadas: ", cidades_visitadas)
   # sol.x = zeros(Int, n, n) # Cria uma matriz de zeros com N linhas e N colunas
    
    

    fechar_ciclo = first(cidades_visitadas)
    insert!(cidades_visitadas, length(cidades_visitadas) + 1, fechar_ciclo)  # Insere a cidade mais próxima após cidade_inicial

    
  
    distancia_total = 0
    for i in 1:(length(cidades_visitadas) - 1)
        distancia_total += data.c[cidades_visitadas[i], cidades_visitadas[i + 1]]
    end 
    sol.x = zeros(Int, n, n) # Cria uma matriz de zeros com N linhas e N colunas
    for i in 1:(length(cidades_visitadas) - 1)
        sol.x[cidades_visitadas[i], cidades_visitadas[i + 1]] = 1
       # println("cidade inicial: ", cidades_visitadas[i], " cidade visitada: ", cidades_visitadas[i + 1]) 
        #println("sol[",cidades_visitadas[i],",",cidades_visitadas[i + 1],"] = ",sol.x[cidades_visitadas[i], cidades_visitadas[i + 1]])
       
    end
    
    sol.route = cidades_visitadas
    sol.primal_bound = distancia_total
    stats.UB_HC = sol.primal_bound
   # println("best_UB: ", stats.UB_HC)
    return
end

function solve_C_GRASP!(data::InstanceData,  data_type::InstanceTypeData,sol::StdFormModelSolution, stats::StatisticsData)

   # println(" Construtive heuristcs GRASP...")

    alfa = 0.25
    n = data_type.DIMENSION
    cidades_visitar = collect(2:n)  # Cria um vetor com N itens que vão de 1 até N
    #sol.x = zeros(Int, n, n) # Cria uma matriz de zeros com N linhas e N colunas
    
    
    
    

    distancia_total = 0
    cidades_visitadas = [1]
    while length(cidades_visitar) > 0

        
        cidades = []#Vector(undef, length(cidades_visitar))
        distancias = []#fill(Inf, length(cidades_visitar))
       
       # cidades_mais_proximas = Vector{Any}(undef, div(length(cidades_visitar), 4))
        #println("cidades mais proximas: ", length(cidades_mais_proximas))
        #distancias_mais_proximas = fill(Inf, length(cidades_mais_proximas))
        
        #println("distancias mais proximas: ", length(distancias_mais_proximas))
        i = cidades_visitadas[end]

        #println("i",i)
       # se cidades_mais_proximas tiver algum elemento vazio distancia_mais_proxima = Inf se não distancia_mais_proxima = menor valor de 

       

        if length(cidades_visitar) != 1

            for j in cidades_visitar #first(cidade):last(cidade)
                if j != i
                    distancia = data.c[i, j]
                    #println("data.c[$i, $j]",data.c[i, j])
                    

                    push!(distancias, distancia)
                    push!(cidades, j)
                    
                end
            end
            #println("distancias: ", distancias)
            #println("cidades: ", cidades)


            distancia_maior = maximum(distancias)
            distancia_menor = minimum(distancias)
            distancia_maxima = distancia_menor + alfa * (distancia_maior - distancia_menor)
            #println("distancia maxima: ", distancia_maxima)
            cidades_mais_proximas = [cidades[i] for i in 1:length(cidades) if distancias[i] <= distancia_maxima]
            distancias_mais_proximas = [distancias[i] for i in 1:length(distancias) if distancias[i] <= distancia_maxima]
            #println("distancias mais proximas: ", distancias_mais_proximas)
            #println("cidades mais proximas: ", cidades_mais_proximas)


            distancia_aleatorio = rand(distancias_mais_proximas)
            pos = findfirst(==(distancia_aleatorio), distancias_mais_proximas)  # Encontra a posição do valor cidade_inicial no vetor cidade
                # println("Posição: ", pos)
            cidade_aleatoria = cidades_mais_proximas[pos]	
           # println("distancia aleatoria: ", distancia_aleatorio)
           # println("cidade aleatoria: ", cidade_aleatoria)
            #insert!(cidades_visitadas, pos , cidade_aleatoria)  # Insere a cidade mais próxima após cidade_inicial
            push!(cidades_visitadas, cidade_aleatoria)    
            distancia_total = distancia_total + distancia_aleatorio
            # println("Distancia total: ", distancia_total)

            if length(cidades_visitadas) > 0
                    pos = findfirst(==(cidade_aleatoria), cidades_visitar)  # Encontra a posição do valor cidade_inicial no vetor cidade
                    # println("Posição: ", pos)

                    deleteat!(cidades_visitar, pos)  # Remove a cidade inicial do vetor cidade
                    # println("Cidades: ", cidade)
            end
       
    
         sol.x[i, cidade_aleatoria] = 1
         #println("sol.x[$i, $cidade_aleatoria]", sol.x[i, cidade_aleatoria])
            # cidade_inicial = cidade_mais_proxima
            
         #println("cidades_visitar", cidades_visitar)
         #println("cidades_visitadas", cidades_visitadas)


        else
            distancia = data.c[cidades_visitadas[end],cidades_visitar[1]] + data.c[cidades_visitar[1],1]
          #  println("entrou aqui")
            push!(cidades_visitadas, cidades_visitar[1])  
            push!(cidades_visitadas, 1)

            
            

            sol.x[cidades_visitadas[end], cidades_visitar[1]] = 1
            sol.x[cidades_visitar[1], 1] = 1
            deleteat!(cidades_visitar, 1)
            # cidade_inicial = cidade_mais_proxima

         #println("cidades_visitar", cidades_visitar)
         #println("cidades_visitadas", cidades_visitadas)

         distancia_total = distancia_total + distancia
         #println("Distancia total: ", distancia_total)
        end

    end  
    
    sol.route = cidades_visitadas
    sol.primal_bound = distancia_total
    stats.UB_HC = sol.primal_bound
   # println("best_UB: ", stats.UB_HC)
    return
end
end # module