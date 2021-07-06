using DataStructures
using DelimitedFiles

"""
This function loads the FUND input parameters stored in datadir into the syntax
needed for the @defmcs macro. This script writes to a text file which can then 
be copied into a @defmcs macro.
"""
function load_mcs_RVs(; txt_out=joinpath(@__DIR__, "mcs_RVs.txt"), datadir=joinpath(@__DIR__, "../../data_for_load_mcs_RVs"), string_labels=false)

    files = readdir(datadir)
    filter!(i -> i != "desktop.ini", files)
    parameters = OrderedDict{Any,Any}(splitext(file)[1] => readdlm(joinpath(datadir, file), ',') for file in files)

    open(txt_out, "w") do f 
        for (k, v) in parameters
            val = loadparameter(v, string_labels)
            if val != nothing
                write(f, "$k = $val\n")
            end
        end 
    end

end

# Helper function for concatenating strings in loadparameter
_concat(x,y) = string(x, ", ", y)

"""
Takes the contents of one parameter file and turns it into the string needed for a definition within a @defmcs macro.
Input 'p' may be a scalar, a 1-D array (defined as two columns in long format), or a 2-D array (defined as three columns in long format).
If p contains no distributional values, then this returns nothing, because it will not be a random variable for mcs.
"""
function loadparameter(p, string_labels::Bool)
    column_count = size(p, 2)
    if column_count == 1
        return format_distribution(p[1,1])
    elseif column_count == 2
        vals = []
        for i in 1:size(p, 1)
            dist = format_distribution(p[i, 2])
            if dist !== nothing 
                string_labels ? push!(vals, "\"$(p[i,1])\" => $dist") : push!(vals, "$(p[i,1]) => $dist")
            end
        end
        return isempty(vals) ? nothing : "[$(reduce(_concat, vals))]"
    elseif column_count == 3
        vals = []
        for i in 1:size(p, 1)
            dist = format_distribution(p[i, 3])
            if dist !== nothing
                string_labels ? push!(vals, "[\"$(p[i,1])\", \"$(p[i,2])\"] => $dist") : push!(vals, "[$(p[i,1]), $(p[i,2])] => $dist")
            end
        end
        return isempty(vals) ? nothing : "[$(reduce(_concat, vals))]"
    else
        error("Unable to parse data files with $column_count columns.")
    end
end

"""
Converts one distribution as defined in FUND's datafiles into a string of a julia-parsable distribution definition.
Example: "~N(0.089;0.1484;min=0)" becomes "Truncated(Normal(0.089,0.1484), 0.0, Inf)"
If v is a scalar or boolean, it returns nothing.
If v starts with '~', but is not ~N(), ~Gamma(), or ~Triangular(), it throws an error. 
"""
function format_distribution(v)
    if isa(v, AbstractString)
        if startswith(v, "~") & endswith(v, ")")
            args_start_index = findnext('(', v, 1)
            dist_name = v[2:args_start_index - 1]
            args = split(v[args_start_index + 1:end - 1], ';')
            fixedargs = filter(i -> !contains(i, "="), args)
            optargs = Dict(split(i, '=')[1] => split(i, '=')[2] for i in filter(i -> contains(i, "="), args))

            if dist_name == "N"
                if length(fixedargs) != 2 error() end
                if length(optargs) > 2 error() end

                basenormal = "Normal($(parse(Float64, fixedargs[1])),$(parse(Float64, fixedargs[2])))"

                if length(optargs) == 0
                    return basenormal
                else
                    min = haskey(optargs, "min") ? parse(Float64, optargs["min"]) : -Inf
                max = haskey(optargs, "max") ? parse(Float64, optargs["max"]) : Inf
                    return "Truncated($basenormal, $min, $max)"
                end
            elseif startswith(v, "~Gamma(")
                if length(fixedargs) != 2 error() end
                if length(optargs) > 2 error() end

                basegamma = "Gamma($(parse(Float64, fixedargs[1])),$(parse(Float64, fixedargs[2])))"

                if length(optargs) == 0
                    return basegamma
                else
                    min = haskey(optargs, "min") ? parse(Float64, optargs["min"]) : -Inf
                    max = haskey(optargs, "max") ? parse(Float64, optargs["max"]) : Inf
                    return "Truncated($basegamma, $min, $max)"
                end
            elseif startswith(v, "~Triangular(")
                triang = "TriangularDist($(parse(Float64, fixedargs[1])), $(parse(Float64, fixedargs[2])), $(parse(Float64, fixedargs[3])))"
                return triang
            else
                error("Unknown distribution")
            end
        end
    end
    return nothing
end
