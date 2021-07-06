using DataFrames, DataFramesMeta, Mimi

include("../src/main.jl")

df_diff = DataFrame([UTF8String, Float64], [:name, :maxdiff], 0)
for c in Mimi.components(m), v in Mimi.variables(m, c)
    filename = "csharp_output/$c.$v.csv"
    if isfile(filename)
        t = typeof(m[c,v])
        if t == Float64
            s = readall(filename)
            csharp_val = parse(Float64, s)
            diff = csharp_val - m[c,v]
            push!(df_diff, (string(v), diff))
        elseif t == Array{Float64,1}
            csharp_val = vec(readdlm(filename))
            jl_val = m[c,v]
            if length(jl_val) != 16
                jl_val = jl_val[1:end - 1]
            end
            diff = csharp_val .- jl_val
            push!(df_diff, (string(v), maxabs(diff)))
        elseif t == Array{Float64,2}
            csharp_val = readcsv(filename)
            jl_val = m[c,v]
            if size(jl_val, 1) != 16
                jl_val = jl_val = m[c,v][1:end - 1,:]
            end
            diff = csharp_val .- jl_val
            push!(df_diff, (string(v), maxabs(diff)))
        else
            println(c, ".", v)
        end
    else
        println(c, ".", v)
    end
end

@where(df_diff, :maxdiff .> 0.00000000001)
