module Formulation

using JuMP
using Gurobi
using Data
using OutputStatistics
using Combinatorics


export create_stdform_model!, solve_stdform_model!


function create_stdform_model!(data::InstanceData,  data_type::InstanceTypeData, model::Model)

    println("Creating model...")

    N = 1:data_type.DIMENSION
    n = Int(data_type.DIMENSION)

    @variable(model, x[i in N, j in N], Bin)
    @variable(model, u[i in N], Int)
 
     #fo
    @objective(model, Min, sum(data.c[i,j] * x[i, j] for i in N, j in N if j != i))

    @constraint(model, [j in N], sum(x[i, j] for i in N  if i != j) == 1)

    @constraint(model, [i in N], sum(x[i, j] for j in N  if j != i) == 1)


    @constraint(model, [i in 2:n, j in 2:n , i !=  j  ], u[i] - u[j] + n * x[i, j] <= n - 1)


    

   #=  @constraint(model, x[1, 14] == 1)
    @constraint(model, x[14, 13] == 1)
    @constraint(model, x[13, 12] == 1)
    @constraint(model, x[12, 7] == 1)
    @constraint(model, x[7, 6] == 1)
    @constraint(model, x[6, 15] == 1)
    @constraint(model, x[15, 5] == 1)
    @constraint(model, x[5, 11] == 1)
    @constraint(model, x[11, 9] == 1)
    @constraint(model, x[9, 10] == 1)
    @constraint(model, x[10, 16] == 1)
    @constraint(model, x[16, 3] == 1)
    @constraint(model, x[3, 2] == 1)
    @constraint(model, x[2, 4] == 1)
    @constraint(model, x[4, 8] == 1)
    @constraint(model, x[8, 1] == 1) =#
    return
end

function solve_stdform_model!(model::Model, data::InstanceData, data_type::InstanceTypeData, sol::StdFormModelSolution, stats::StatisticsData) 

    N = 1:data_type.DIMENSION
    optimize!(model)

    sol.status = termination_status(model)
    sol.primal_bound = 1e8
    sol.dual_bound = -1e8

    fill!(sol.x, 0.0)

    stats.sol_status = sol.status

    if has_values(model)

        sol.primal_bound = objective_value(model)
        sol.dual_bound = objective_bound(model)

        stats.total_time = solve_time(model)
        stats.best_LB = sol.dual_bound
        stats.best_UB = sol.primal_bound
        stats.gap = 100 * ((stats.best_UB - stats.best_LB) / stats.best_UB)
        
        for i in N, j in N
            
            if value(model[:x][i,j]) * data.c[i,j]  > 0.0001
                sol.x[i, j] = value(model[:x][i, j])
            end

        end

    end

end

end # module