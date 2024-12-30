module Disturbance

using JuMP
using Data
using OutputStatistics
using Random


export DoubleBridge!

function DoubleBridge!(data::InstanceData, data_type::InstanceTypeData, sol::StdFormModelSolution, stats::StatisticsData)

  #println(" pertubação DoubleBridge ...")


  n = data_type.DIMENSION
  # é quantas cidades tem no total, o tamanho da rota é n+1

  rota_inicial = sol.route
  #println("rota_inicial: ", rota_inicial)

  # Gera um índice inicial aleatório entre 1 e length(rota_inicial) - 1
  
  indice_inicial1 = rand(2:round(Int, (n/3)-1)) 
  #println("indice_inicial: ", indice_inicial1)
  segmento1 = rota_inicial[indice_inicial1:indice_inicial1 + 1]
  
  #println("Segmento 1: ", segmento1)
  tamanho_segmento2 = round(Int, (n/3)-1)
  #println("tamanho_segmento2: ", tamanho_segmento2)
  
  indice_inicial2 = 0

  indice_inicial2 = rand(round(Int, (n/3))+2:n-tamanho_segmento2)
  #println("Índice inicial segmento2: ", indice_inicial2)
  #= while true
      
      # indice_inicial2 = rand(round(Int, (n/3)-1)+2:n-1)
      if indice_inicial2 < indice_inicial1 - tamanho_segmento2 -1 || indice_inicial2 > indice_inicial1 + 1
        #   println("Índice inicial segmento2: ", indice_inicial2)
          break
      end

  end =#

  segmento2 = rota_inicial[indice_inicial2:indice_inicial2 + tamanho_segmento2 - 1]
  
  
  #println("Segmento 2: ", segmento2)


  delta = 0
  delta -= data.c[rota_inicial[indice_inicial1-1],rota_inicial[indice_inicial1]]
  # println("data.c[$(rota_inicial[indice_inicial1-1]),$(rota_inicial[indice_inicial1])]: ")
  delta -= data.c[rota_inicial[indice_inicial1+1],rota_inicial[indice_inicial1+2]]
  # println("data.c[$(rota_inicial[indice_inicial1+1]),$(rota_inicial[indice_inicial1+2])]: ")
  delta -= data.c[rota_inicial[indice_inicial2-1],rota_inicial[indice_inicial2]]
  # println("data.c[$(rota_inicial[indice_inicial2-1]),$(rota_inicial[indice_inicial2])]: ")
  delta -= data.c[rota_inicial[indice_inicial2+tamanho_segmento2-1],rota_inicial[indice_inicial2+tamanho_segmento2]]
  #println("data.c[$(rota_inicial[indice_inicial2+tamanho_segmento2-1]),$(rota_inicial[indice_inicial2+tamanho_segmento2])]: ")
  delta += data.c[rota_inicial[indice_inicial1-1],rota_inicial[indice_inicial2]]
  #println("data.c[$(rota_inicial[indice_inicial1-1]),$(rota_inicial[indice_inicial2])]: ")


  delta += data.c[rota_inicial[indice_inicial2+tamanho_segmento2-1],rota_inicial[indice_inicial1+2]]
  #println("data.c[$(rota_inicial[indice_inicial2+tamanho_segmento2-1]),$(rota_inicial[indice_inicial1+2])]: ")


  delta += data.c[rota_inicial[indice_inicial2-1],rota_inicial[indice_inicial1]]
  #println("data.c[$(rota_inicial[indice_inicial2-1]),$(rota_inicial[indice_inicial1])]: ")
  delta += data.c[rota_inicial[indice_inicial1+1],rota_inicial[indice_inicial2+tamanho_segmento2]]
  #println("data.c[$(rota_inicial[indice_inicial1+1]),$(rota_inicial[indice_inicial2+tamanho_segmento2])]: ")

  # Remove segmento 1
  pos1 = findfirst(==(rota_inicial[indice_inicial1]), rota_inicial)  # Encontra a posição do valor cidade_inicial no vetor cidade
  elemento = rota_inicial[indice_inicial1-1]
  #println("elemento: ", elemento)
  deleteat!(rota_inicial, pos1)  # Remove a cidade inicial do vetor cidade
  deleteat!(rota_inicial, pos1)  # Remove a cidade inicial do vetor cidade
  insert!(rota_inicial, pos1, elemento)  # Insere a cidade mais próxima após cidade_inicial
  #println("rota_inicial após remoção do segmento 1: ", rota_inicial)
  
  
  
  pos2 = findfirst(==(rota_inicial[indice_inicial2]), rota_inicial)
  insert!(rota_inicial, pos2, segmento1[1])  # Insere a cidade mais próxima após cidade_inicial
  insert!(rota_inicial, pos2 + 1, segmento1[2])  # Insere a cidade mais próxima após cidade_inicial
  # println("rota_inicial após readicionar do segmento 1: ", rota_inicial)


  # Remove segmento 2
  for i in segmento2
    pos = findfirst(==(i), rota_inicial)
    deleteat!(rota_inicial, pos)
      
  end
  #  println("rota_inicial após remoção do segmento 2: ", rota_inicial)


  pos3 = findfirst(==(elemento), rota_inicial)

  splice!(rota_inicial, pos3+1, segmento2)

  #println("rota_inicial após remoção do segmento 2: ", rota_inicial)


  #println("rota_inicial após pertubação DoubleBridge: ", rota_inicial)

  sol.route = rota_inicial

  sol.primal_bound = sol.primal_bound + delta 
      
  stats.UB_LS = sol.primal_bound
      

  #println("best_UB_LS: ", stats.UB_LS)



       
       
   

    
    


end

end # module