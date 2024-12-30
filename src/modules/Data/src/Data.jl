module Data

using DelimitedFiles

struct InstanceTypeData
    instName::String
    DIMENSION:: Int
    COORD_Type::String
    EDGE_WEIGHT_FORMAT::String
    Aux:: Int
end

struct InstanceData
    c::Array{Float64}  # Cost
end

export InstanceTypeData, readTypeData, InstanceData, readDataXY, readDataEXPLICIT

function readTypeData(instanceFile::String)

    file = open(instanceFile)
    fileText = read(file, String)
    tokens = split(fileText) 
    #tokens will have all the tokens of the input 
    #file in a single vector. We will get the input token by token

 
   
    
    aux = 1
    i = 1
    name  = ""
    dim = ""
    coord =""
    weights = "sem"
    while aux == 1

        if tokens[i] == "NAME:"
            name = tokens[i+1]   
        end

        if tokens[i] == "DIMENSION:"
            dim = parse(Int,tokens[i+1])    
        end

        if tokens[i] == "EDGE_WEIGHT_TYPE:"
            coord = tokens[i+1]
        end
        if coord == "EXPLICIT"
            if tokens[i] == "EDGE_WEIGHT_FORMAT:"
                weights = tokens[i+1]                
            end
        end
        if tokens[i] == "NODE_COORD_SECTION" || tokens[i] == "EDGE_WEIGHT_SECTION"
            aux = i
            
        end
        i += 1
    end
    
  

    # Print instance data
    print("Instance TYPE: ", 
            "\n",name ,
            "\nDIMENSION: ",dim,
            "\nEDGE_WEIGHT_TYPE: ", coord,
            "\nEDGE_WEIGHT_FORMAT:",weights,
            "\ni= ",aux,                       #we will start reading the data from here
            "\ntamnho ",length(tokens))         #length tokens

           

        datatype = InstanceTypeData(name,dim,coord,weights,aux)

    return datatype

    
  

end


function readDataXY(instanceFile::String, datatype::InstanceTypeData)

    file = open(instanceFile)

    fileText = read(file, String)
    tokens = split(fileText) 
    

    i = datatype.Aux
    j = 1
   
    dim = datatype.DIMENSION

    
    x = Array{Float32}(undef,dim)
    y = Array{Float32}(undef,dim)

    if datatype.COORD_Type == "CEIL_2D"
        while i < length(tokens)-3     
            
            x[j] = ceil(parse(Int64,tokens[i+2]) )
        
            y[j] = ceil(parse(Int64,tokens[i+3]) )
            j += 1
            i += 3   

        end
    else
        while i < length(tokens)-3   
            x[j] = parse(Float32,tokens[i+2]) 
            
            y[j] = parse(Float32,tokens[i+3]) 
            j += 1
            i += 3 
        end  
    end
 
    if datatype.COORD_Type == "ATT"
        println("ATT")
        xd = Array{Float32}(undef,dim,dim)
        yd = Array{Float32}(undef,dim,dim)
        r = Array{Float64}(undef, dim, dim)
        c = Array{Float64}(undef, dim, dim)

        for i = 1:dim
            for j = 1:dim
                if i <= j
                    xd[i,j] = (x[i] - x[j])^2
                    yd[i,j] = (y[i] - y[j])^2
                    r[i,j] = r[j,i] = ceil(sqrt((xd[i,j] +yd[i,j] )/10.0))
                end
             

            end
        end 
  
    elseif datatype.COORD_Type == "GEO"
       
      
       

       
        println("GEO")
        PI = 3.141592
       
       

        degx = Array{Float64}(undef,dim)
        degy = Array{Float64}(undef,dim)
        minx = Array{Float64}(undef,dim)
        miny = Array{Float64}(undef,dim)
        latitude = Array{Float64}(undef,dim)
        longitude = Array{Float64}(undef,dim)
        for i = 1:dim
          
            
            degx[i] = trunc(x[i])
            degy[i] = trunc(y[i])
            minx[i] = x[i] - degx[i]
            miny[i] = y[i] - degy[i]
            latitude[i] = (PI * (degx[i] + (5 * minx[i] )/ 3))/180
            longitude[i] = (PI * (degy[i] + (5 * miny[i]) / 3))/180
    
        end
        
        
        RRR = 6378.388;
        q1 = Array{Float64}(undef,dim,dim)
        q2 = Array{Float64}(undef,dim,dim)
        q3 = Array{Float64}(undef,dim,dim)
        for i = 1:dim
            for j = 1:i
            
                q1[i,j] = cos( longitude[i] - longitude[j] )
                q2[i,j]  = cos( latitude[i] - latitude[j] )
                q3[i,j]  = cos( latitude[i] + latitude[j] )
            end
        end
        
        c = zeros(Float64, dim, dim)

        for i = 1:dim
            for j = 1:i
               
                   #      RRR.*acos(0.5.*((1.0.+q1       )  .  *q2.    -(1.0.       -q1).*q3       )).+1.0
                c[i,j] = trunc(RRR * acos(0.5 *((1.0 + q1[i,j] ) * q2[i,j]  - (1.0 - q1[i,j] )* q3[i,j]  )) + 1.0)


            end
        end 


        for j = 1:dim
            for i = j+1:dim
               
                
                c[j,i] =  c[i,j] 


            end
        end 
        
    elseif datatype.COORD_Type == "EUC_2D" 

        xd = Array{Float32}(undef,dim,dim)
        yd = Array{Float32}(undef,dim,dim)
        c = Array{Float64}(undef, dim, dim)

        for i = 1:dim
            for j = 1:i
                
                xd[i,j] = x[i] - x[j]
                yd[i,j] = y[i] - y[j]
                c[i,j] = c[i,j] = round(sqrt((xd[i,j]^2 + yd[i,j]^2)),RoundNearestTiesUp)
                
            end
        end

        for j = 1:dim
            for i = j+1:dim
               
                
                c[j,i] =  c[i,j] 


            end
        end 
    elseif datatype.COORD_Type == "CEIL_2D" 

        xd = Array{Float32}(undef,dim,dim)
        yd = Array{Float32}(undef,dim,dim)
        c = Array{Float64}(undef, dim, dim)

        for i = 1:dim
            for j = 1:i
                
                xd[i,j] = x[i] - x[j]
                yd[i,j] = y[i] - y[j]
                c[i,j] = c[j,i] = ceil(sqrt((xd[i,j]^2 + yd[i,j]^2)))
                
            end
        end
        for j = 1:dim
            for i = j+1:dim
               
                
                c[j,i] =  c[i,j] 


            end
        end 

    end
    
    
 
    instance = InstanceData(c)

    return instance

end

function readDataEXPLICIT(instanceFile::String, datatype::InstanceTypeData)
   
    file = open(instanceFile)
    fileText = read(file, String)
    tokens = split(fileText) 
    

    i = datatype.Aux +1
    j = 1
    dim = datatype.DIMENSION
    #c = Array{Float64}(undef, dim, dim)
    c = zeros(Float32,dim,dim)

    if datatype.EDGE_WEIGHT_FORMAT == "LOWER_DIAG_ROW"
        for k = 1:dim
            for j = 1:k

                c[k,j] = parse(Float32,tokens[i]) 
            
                i += 1
            end 

        end
        for k = 1:dim
            for j = k:dim

                c[k,j] = c[j,k]
            
               
            end 

        end
    elseif datatype.EDGE_WEIGHT_FORMAT == "UPPER_ROW"
        for k = 1:dim
            for j = k:dim
                if k != j
                    c[k,j] = parse(Float32,tokens[i]) 
                    #println("custo C[$k,$j] = ", c[k,j])
                    i += 1
                end
            end 

            for k = 1:dim
                for j = 1:k
    
                    c[k,j] = c[j,k]
                
                   
                end 
    
            end

        end
    elseif datatype.EDGE_WEIGHT_FORMAT == "FULL_MATRIX"
        for k = 1:dim
            for j = 1:dim
               
                c[k,j] = parse(Float32,tokens[i]) 
                #println("custo C[$k,$j] = ", c[k,j])
                i += 1
               
            end 
        end

    elseif datatype.EDGE_WEIGHT_FORMAT == "UPPER_DIAG_ROW"
        for k = 1:dim
            for j = k:dim
               
                c[k,j] = parse(Float32,tokens[i]) 
                #println("custo C[$k,$j] = ", c[k,j])
                i += 1
               
            end 
        end

        for k = 1:dim
            for j = 1:k

                c[k,j] = c[j,k]
            
               
            end 

        end

    end

    

    instance = InstanceData(c)

    return instance

end
end # module